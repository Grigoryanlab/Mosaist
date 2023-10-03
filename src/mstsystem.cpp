#include "mstsystem.h"

using namespace MST;

string MstSys::pathBase(const string& fn) {
  if (fn.find_last_of(".") == string::npos) return fn;
  else return fn.substr(0, fn.find_last_of("."));
}

string MstSys::splitPath(const string& path, int outToken, string* dirPathPtr, string* fileNamePtr, string* extensionPtr) {
  string dirPath, fileName, extension;
  int pos = path.rfind("/");
  if (pos == string::npos) {
    fileName = path;
    dirPath = "./";
  } else if (pos == 0) {
    fileName = path.substr(1);
    dirPath = "/";
  } else {
    fileName = path.substr(pos+1);
    dirPath = path.substr(0, pos);
  }
  pos = fileName.rfind(".");
  if (pos == string::npos) {
    extension = "";
  } else {
    extension = fileName.substr(pos+1);
    fileName = fileName.substr(0, pos);
  }
  if (dirPathPtr != NULL) *dirPathPtr = dirPath;
  if (fileNamePtr != NULL) *fileNamePtr = fileName;
  if (extensionPtr != NULL) *extensionPtr = extension;
  switch (outToken) {
    case 0:
      return dirPath;
    case 1:
      return fileName;
    case 2:
      return extension;
    default:
      MstUtils::error("unrecognized output token type specified '" + MstUtils::toString(outToken) + "'", " MstSys::splitPath");
  }
  return ""; // just to make the compiler happy, this is never reached
}

bool MstSys::fileExists(const char* filename) {
  struct stat buffer ;
  if (stat( filename, &buffer) == 0) return true;
  return false;
}

long MstSys::fileSize(const char* filename) {
  struct stat stat_buf;
  int rc = stat(filename, &stat_buf);
  return rc == 0 ? stat_buf.st_size : -1;
}

bool MstSys::isDir(const char *filename) {
  struct stat buffer ;
  if (stat( filename, &buffer) < 0) return false;
  return (buffer.st_mode & S_IFDIR);
}

int MstSys::csystem(const string& cmd, bool checkError, int success, const string& from) {
  int ret = system(cmd.c_str());
  string head = from.empty() ? "" : from + " -> ";
  if ((checkError) && (ret != success)) {
    MstUtils::error("system command '" + cmd + "' failed", head + "MstUtils::csystem", ret);
  }
  return ret;
}

void MstSys::cmkdir(const string& dirPath, bool makeParents) {
  if (MstSys::isDir(dirPath)) return; // return if the path is already a dir
  int ret = MstSys::csystem("mkdir " + (makeParents ? (string) "-p " : (string) "") + dirPath, false);
  MstUtils::assertCond(ret == 0, "failed to make directory '" + dirPath + "'");
}

void MstSys::crmdir(const string& dirPath, bool recursive) {
  int ret = MstSys::csystem((recursive ? (string) "rmdir " : (string) "rm -r ") + dirPath, false);
  MstUtils::assertCond(ret == 0, "failed to remove directory '" + dirPath + "'" + (recursive ? " recursively" : ""));
}

void MstSys::crm(const string& filePath) {
  int ret = MstSys::csystem("rm " + filePath, false);
  MstUtils::assertCond(ret == 0, "failed to remove file '" + filePath + "'");
}

string MstSys::getMachineName() {
  int n = 1024;
  char hostname[n];
  hostname[n-1] = '\0';
  gethostname(hostname, n-1);
  return string(hostname);
}

string MstSys::getUserName() {
  int n = 1024;
  char username[n];
  username[n-1] = '\0';
  getlogin_r(username, n-1);
  return string(username);
}

bool MstSys::getNetLock(const string& tag, bool shared, const string& linuxHost) {
  int timeout = 60*5, waitTime = 10; // in seconds
  string eQ = "\'\"\'\"\'";
  string token = getMachineName() + "-" + MstUtils::toString((int) getpid());
  string timeoutStr = MstUtils::toString(timeout);
  string waitStr = MstUtils::toString(waitTime);
  string base = "/tmp/.mst-" + tag;
  string lockFile = base + ".lockfile";
  string busyFileX = base + ".busyfile.X";
  string busyFileS = base + ".busyfile.S";
  string busyFiles = busyFileX + " " + busyFileS;
  string notLockedCond = shared ? "[ ! -f " + busyFileX + " ]" : "[ ! -f " + busyFileX + " ] && [ ! -f " + busyFileS + " ] ";
  string doLock = shared ? "echo " + token + " >> " + busyFileS : "echo " + token + " > " + busyFileX;
  string timeCheck = "if ( [ ! -f " + busyFileX + " ] || [ $((`date +%s` - `stat -c %Y " + busyFileX + "`)) -gt " + timeoutStr + " ] ) && "
                        "( [ ! -f " + busyFileS + " ] || [ $((`date +%s` - `stat -c %Y " + busyFileS + "`)) -gt " + timeoutStr + " ] ); then "
                          "rm -f " + busyFiles + "; " + doLock + "; exit 3; else exit 2; fi";
  string touchCmd = "if " + notLockedCond + "; then " + doLock + " ; exit 0; else " + timeCheck + "; fi";
  string remoteCmd = "flock -w " + waitStr + " " + lockFile + " -c " + eQ + touchCmd + eQ;
  string cmd = "ssh " + linuxHost + " '" + remoteCmd + "'";
  for (int i = 0; i < int(1.0*timeout/waitTime + 1); i++) {
    int ret = system(cmd.c_str());
    ret = WEXITSTATUS(ret);    // this gets the exit status of the command
    if (ret == 0) return true; // got lock fine
    if (ret == 3) return true; // got lock, but got impatient and had to overwrite the old one
    if (ret == 2) sleep(waitTime); // could not get lock, somebody else's lock is not yet released
    // (ret == 1) means getting a lock failed; not clear how to recover from that
  }
  return false;
}

bool MstSys::releaseNetLock(const string& tag, const string& linuxHost) {
  int timeout = 60*5, waitTime = 10; // in seconds
  string eQ = "\'\"\'\"\'";
  string token = getMachineName() + "-" + MstUtils::toString((int) getpid());
  string timeoutStr = MstUtils::toString(timeout);
  string waitStr = MstUtils::toString(waitTime);
  string base = "/tmp/.mst-" + tag;
  string lockFile = base + ".lockfile";
  string busyFileX = base + ".busyfile.X";
  string busyFileS = base + ".busyfile.S";
  string sed = "sed -i " + eQ + "/^" + token +"$/d" + eQ + " ";
  string removeTokenX = "if [ -f " + busyFileX + " ]; then " + sed + busyFileX + "; fi; ";
  string removeTokenS = "if [ -f " + busyFileS + " ]; then " + sed + busyFileS + "; fi; ";
  string removeToken = removeTokenX + removeTokenS; // there should never be BOTH an exclusive and a shared lock
  string cleanX = "if [ ! -s " + busyFileX + " ]; then rm -f " + busyFileX + "; fi; ";
  string cleanS = "if [ ! -s " + busyFileS + " ]; then rm -f " + busyFileS + "; fi; ";
  string clean = cleanX + cleanS;
  string remoteCmd = "flock -w " + waitStr + " " + lockFile + " -c \"" + removeToken + clean + "\"; exit 0;";
  string cmd = "ssh " + linuxHost + " '" + remoteCmd + "'";
  for (int i = 0; i < int(1.0*timeout/waitTime + 1); i++) {
    int ret = system(cmd.c_str());
    ret = WEXITSTATUS(ret);    // this gets the exit status of the command
    if (ret == 0) return true; // got lock fine
    // (ret == 1) means getting a lock failed; not clear how to recover from that
  }
  return false;
}

int MstSys::memUsage() {
  string pid = MstUtils::toString((int) getpid());
  string tmpFile = "/tmp/mu.out." + pid;
  MstSys::csystem("ps -p " + pid + " -o rss | tail -1 > " + tmpFile);
  vector<string> lines = MstUtils::fileToArray(tmpFile);
  if ((lines.size() != 1) || (!MstUtils::isInt(lines[0]))) MstUtils::error("could not get memory usage, output in " + tmpFile);
  MstSys::crm(tmpFile);
  return MstUtils::toInt(lines[0]);
}
