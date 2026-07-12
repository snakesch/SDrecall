#include <time.h>

/* Keep timestamp-bearing tool headers reproducible across sequential runs. */
time_t time(time_t *result)
{
    const time_t fixed_time = 1783814400;

    if (result != NULL) {
        *result = fixed_time;
    }
    return fixed_time;
}
