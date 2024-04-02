#include <sys/stat.h>
#include <sys/types.h>

int c_mkdir755(const char *path) {
    return mkdir(path, 0755);
}

void c_mkdir755noreturn(const char *path) {
    mkdir(path, 0755);
}
