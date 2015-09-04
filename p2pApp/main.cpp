
#include <iocsh.h>

#define epicsExportSharedSymbols
#include "pva2pva.h"

int main(int argc, char *argv[])
{
    registerGWClientIocsh();
    registerGWServerIocsh();
    if(argc>1)
        iocsh(argv[1]);
    int ret = iocsh(NULL);
    gwServerShutdown();
    gwClientShutdown();
    return ret;
}
