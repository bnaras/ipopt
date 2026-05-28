#include <IpStdCInterface.h>

int main(void) {
  int major = 0, minor = 0, release = 0;
  GetIpoptVersion(&major, &minor, &release);
  return major < 3;
}
