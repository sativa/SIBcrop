#ifndef SIB_NC_UTIL_H
#define SIB_NC_UTIL_H

#define CHECK(expr) call nc_check( expr, __FILE__, __LINE__)

#define ENSURE_VAR(ncid,varname,varid) \
    call nc_ensure_var(ncid, varname, varid, __FILE__, __LINE__)

#endif
