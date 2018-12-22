/*
 * userSettings.h
 *
 *  Created on: Jan 20, 2014
 *      Author: Claas Grohnfeldt
 */

#ifndef USERSETTINGS_H_
#define USERSETTINGS_H_


#include "includes.h"

// project headers
#include "dataIO.h"
#include "paths.h"

void getUserSettings(SpEOPaths *paths, SpEODataIOSetting *dSetting, SpEOFusionSetting *fSetting, SpEOOutputSetting *oSetting, SpEOSolverSetting *sSetting, SpEOParallelSetting *pSetting, int argc, char **argv);

#endif /* USERSETTINGS_H_ */
