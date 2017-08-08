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
  MstUtils::assert(ret == 0, "failed to make directory '" + dirPath + "'");
}

void MstSys::crmdir(const string& dirPath, bool recursive) {
  int ret = MstSys::csystem((recursive ? (string) "rmdir " : (string) "rm -r ") + dirPath, false);
  MstUtils::assert(ret == 0, "failed to remove directory '" + dirPath + "'" + (recursive ? " recursively" : ""));
}

void MstSys::crm(const string& filePath) {
  int ret = MstSys::csystem("rm " + filePath, false);
  MstUtils::assert(ret == 0, "failed to remove file '" + filePath + "'");
}
