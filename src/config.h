#ifndef CONFIG_H
#define CONFIG_H

/** \file config.h Routines for manipulating the configuration YML file. */

/** @{ */
/** Check configuration.yml for documentation. */
extern char *HP_CHAIN;
extern int EPS_HH;
extern int EPS_HP;
extern int EPS_HB;
extern int EPS_PP;
extern int EPS_PB;
extern int EPS_BB;
extern int PENALTY_VALUE;
extern int N_CYCLES;
extern int COLONY_SIZE;
extern double FORAGER_RATIO;
extern int IDLE_LIMIT;
extern int N_HIVES;
extern int RANDOM_SEED;
/** @} */

/** Initializes configuration based on the configuration file. */
void initialize_configuration();

#endif // CONFIG_H
