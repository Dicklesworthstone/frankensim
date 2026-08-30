#ifndef FRANKENSIM_APPLE_H
#define FRANKENSIM_APPLE_H

#include <stdint.h>

uint32_t frankensim_apple_schema_version(void);
uint64_t frankensim_apple_run(uint32_t experiment_id, double quality, uint32_t seed);
uint64_t frankensim_apple_result_len(void);
double frankensim_apple_result_value(uint64_t index);
int32_t frankensim_apple_last_error(void);

#endif
