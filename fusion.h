#pragma once

#include <ap_fixed.h>
#include <hls_math.h>
#include <hls_vector.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include "hls_stream.h"
#include "ap_int.h"
#include <math.h>

#define STATE_SIZE 22
#define IMU_SAMPLES 160

typedef float fp;

//void fusion
//(
//		fp state[22],
//		fp P[22][22],
//		fp accel_block[160][3],
//		fp gyro_block[160][3],
//		fp gps_pos[3],
//		fp gps_vel[3],
//		fp orient_block[160][4],
//		fp pos_block[160][3]
//);
void fusion
(
		fp state[22],
		fp P[22][22],
		fp accel_block[160][3],
		fp gyro_block[160][3],
		fp mag_last[3],
		fp gps_pos[3],
		fp gps_vel[3],
		fp orient_block[160][4],
		fp pos_block[160][3]
);
