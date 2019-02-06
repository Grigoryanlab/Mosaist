#ifndef _MSTSYSTEM_H
#define _MSTSYSTEM_H

#include "msttypes.h"
#include <unistd.h>
#include <limits.h>


/* A bunch of static routines to perform various OS-related operations */
class MstSys {
  public:
    static string pathBase(const string& fn); // gets the base name of the path (removes the extension)

    /* Returns the directory part of the path, the file name, or the file extension,
     * depending on whether outToken is 0, 1, or 2, respectively. */
    static string splitPath(const string& path, int outToken, string* dirPathPtr = NULL, string* fileNamePtr = NULL, string* extensionPtr = NULL);

    static bool fileExists(const char *filename);
    static bool fileExists(const string filename) { return fileExists(filename.c_str()); }
    static long fileSize(const char* filename);
    static long fileSize(const string filename) { return fileSize(filename.c_str()); }
    static bool isDir(const char *filename);
    static bool isDir(const string& filename) { return isDir(filename.c_str()); }
    static int csystem(const string& cmd, bool checkError = true, int success = 0, const string& from = "");
    static void cmkdir(const string& dirPath, bool makeParents = false);
    static void crmdir(const string& dirPath, bool recursive = false);
    static void crm(const string& filePath);
    static string getMachineName();
    static string getUserName();
    static bool getNetLock(const string& tag, bool shared = false, const string& linuxHost = "anthill.cs.dartmouth.edu");
    static bool releaseNetLock(const string& tag, const string& linuxHost = "anthill.cs.dartmouth.edu");
    static int memUsage(); // total memory used in Kbytes

  private:
};

#endif
