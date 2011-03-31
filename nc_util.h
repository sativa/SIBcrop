#ifndef SIB_NC_UTIL_H
#define SIB_NC_UTIL_H

#define ENSURE_VAR(ncid,varname,varid) \
    call ensure_var(ncid, varname, varid, __FILE__, __LINE__)

#endif
