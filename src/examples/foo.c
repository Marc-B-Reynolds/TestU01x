
#include <stdint.h>

main(int argc, char** argv)
{
  uint64_t y = (uint64_t)(-1);

  y >>= 12;
  double r = y * 0x1p-52;
  //y >>= 11;
  //double r = y * 0x1p-53;

  //y >>= 1;
  //double r = y * (2.0/(uint64_t)(-1<<10));

  printf("%f %e\n", r, 1-r);
}
