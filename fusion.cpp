#include "fusion.h"

//===========================================================================//
static void quat_normalize(fp q[4])
{
#pragma HLS INLINE
	fp qq0 = q[0]*q[0];
	fp qq1 = q[1]*q[1];
	fp qq2 = q[2]*q[2];
	fp qq3 = q[3]*q[3];
	fp qq = qq0 + qq1 + qq2 + qq3;
	fp n = sqrtf(qq);
	q[0]/=n;
	q[1]/=n;
	q[2]/=n;
	q[3]/=n;
}
//===========================================================================//
static void quat_to_rot(fp q[4], fp R[3][3]) {
    #pragma HLS INLINE
    fp q0 = q[0], q1 = q[1], q2 = q[2], q3 = q[3];
    fp q00 = q0*q0, q11 = q1*q1, q22 = q2*q2, q33 = q3*q3;
    fp q01 = q0*q1, q02 = q0*q2, q03 = q0*q3;
    fp q12 = q1*q2, q13 = q1*q3, q23 = q2*q3;

    R[0][0] = q00 + q11 - q22 - q33;
    R[0][1] = 2*(q12 - q03);
    R[0][2] = 2*(q13 + q02);

    R[1][0] = 2*(q12 + q03);
    R[1][1] = q00 - q11 + q22 - q33;
    R[1][2] = 2*(q23 - q01);

    R[2][0] = 2*(q13 - q02);
    R[2][1] = 2*(q23 + q01);
    R[2][2] = q00 - q11 - q22 + q33;
}

//===========================================================================//
static void inv3x3(fp A[3][3], fp Ainv[3][3])
{
#pragma HLS INLINE
	fp det =
			A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])
			- A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
			+ A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);

	fp invdet = 1.0f/det;

	Ainv[0][0] =  (A[1][1]*A[2][2]-A[1][2]*A[2][1])*invdet;
	Ainv[0][1] = -(A[0][1]*A[2][2]-A[0][2]*A[2][1])*invdet;
	Ainv[0][2] =  (A[0][1]*A[1][2]-A[0][2]*A[1][1])*invdet;

	Ainv[1][0] = -(A[1][0]*A[2][2]-A[1][2]*A[2][0])*invdet;
	Ainv[1][1] =  (A[0][0]*A[2][2]-A[0][2]*A[2][0])*invdet;
	Ainv[1][2] = -(A[0][0]*A[1][2]-A[0][2]*A[1][0])*invdet;

	Ainv[2][0] =  (A[1][0]*A[2][1]-A[1][1]*A[2][0])*invdet;
	Ainv[2][1] = -(A[0][0]*A[2][1]-A[0][1]*A[2][0])*invdet;
	Ainv[2][2] =  (A[0][0]*A[1][1]-A[0][1]*A[1][0])*invdet;
}

//===========================================================================//
void imu_predict
(
		int k,
		fp state[22],
		fp P[22][22],
		fp accel_block[160][3],
		fp gyro_block[160][3],
		fp orient_block[160][4],
		fp pos_block[160][3]
)
{
#pragma HLS ARRAY_PARTITION variable=P complete dim = 2

	static fp Ts = 0.00625;
	static fp g_ned[3] = {0, 0,9.80665};
	fp q[4] = {state[0],state[1],state[2],state[3]};
	fp pos[3] = {state[4],state[5],state[6]};
	fp vel[3] = {state[7],state[8],state[9]};
	fp b_g[3] = {state[10],state[11],state[12]};
	fp b_a[3] = {state[13],state[14],state[15]};

	fp omega[3];
	loop7 : for(int i=0;i<3;i++)omega[i]=gyro_block[k][i]-b_g[i];

	fp dtheta[3];
	loop8 : for(int i=0;i<3;i++)dtheta[i]=omega[i]*Ts;

	fp half_theta[3];

	loop9 : for(int i=0;i<3;i++)
		half_theta[i]=0.5f*dtheta[i];

	fp half_norm2=0;
	loop10 : for(int i=0;i<3;i++)
//#pragma HLS PIPELINE II=3
#pragma HLS UNROLL
		half_norm2+=half_theta[i]*half_theta[i];

	fp dq_w=1-0.5f*half_norm2;

	fp q_new[4];

	q_new[0]= dq_w*q[0] - half_theta[0]*q[1] - half_theta[1]*q[2] - half_theta[2]*q[3];
	q_new[1]= dq_w*q[1] + half_theta[0]*q[0] + half_theta[2]*q[2] - half_theta[1]*q[3];
	q_new[2]= dq_w*q[2] - half_theta[2]*q[1] + half_theta[1]*q[0] + half_theta[0]*q[3];
	q_new[3]= dq_w*q[3] + half_theta[2]*q[0] - half_theta[0]*q[2] + half_theta[1]*q[1];

	quat_normalize(q_new);
	fp Rb2n[3][3];
	quat_to_rot(q_new,Rb2n);
	loop11 : for(int i=0;i<4;i++)
		orient_block[k][i]=q_new[i];

	fp f_body[3];
	loop12 : for(int i=0;i<3;i++)
		f_body[i]=accel_block[k][i]-b_a[i];

	fp f_ned[3] = {0,0,0};
	loop13 : for(int i=0;i<3;i++)
		loop14 : for(int j=0;j<3;j++)
#pragma HLS PIPELINE II=3
			f_ned[i] = f_ned[i] + Rb2n[i][j]*f_body[j];

	fp a_ned[3];
	loop15 : for(int i=0;i<3;i++)
		a_ned[i] = f_ned[i] + g_ned[i];

	fp vel_new[3];
	loop16 : for(int i=0;i<3;i++)
		vel_new[i]= vel[i] + a_ned[i]*Ts;
	fp vel1[3];
	loop17 : for(int i=0;i<3;i++)
#pragma HLS UNROLL
		vel1[i] = vel_new[i] - vel[i];

	fp vel2[3];
	loop18 : for(int i=0;i<3;i++)
		vel2[i] = vel[i] + 0.5*vel1[i];

	loop19 : for(int i=0;i<3;i++)
//#pragma HLS PIPELINE II=3
#pragma HLS UNROLL
		pos[i] = pos[i] + vel2[i] * Ts;

	loop20 : for(int i=0;i<3;i++)
		pos_block[k][i]=pos[i];

	state[0] =q_new[0];
	state[1] =q_new[1];
	state[2] =q_new[2];
	state[3] =q_new[3];
	state[4] =pos[0];
	state[5] =pos[1];
	state[6] =pos[2];
	state[7] =vel_new[0];
	state[8] =vel_new[1];
	state[9] =vel_new[2];

	fp diag_increase[22] = {1e-06,1e-06,1e-06,1e-06,0.02,0.02,0.02,0.05,0.05,0.05,0.0001,0.0001,0.0001,0.005,0.005,0.005,1e-06,1e-06,1e-06,1e-06,1e-06,1e-06};
#pragma HLS ARRAY_PARTITION variable=diag_increase complete dim = 1
	for(int i=0;i<22;i++)
#pragma HLS PIPELINE II=1
		P[i][i] = P[i][i] + diag_increase[i];

}
//===========================================================================//
static void chol6(fp A[6][6], fp L[6][6]) {
    #pragma HLS INLINE
    for(int i=0;i<6;i++) {
        #pragma HLS PIPELINE II=1
        for(int j=0;j<=i;j++) {
            fp sum = 0;
            for(int k=0;k<j;k++) {
                #pragma HLS UNROLL
                sum += L[i][k]*L[j][k];
            }
            L[i][j] = (i == j) ? sqrtf(A[i][i]-sum) : (A[i][j]-sum)/L[j][j];
        }
    }
}

//===========================================================================//
static void chol_solve6(fp S[6][6], fp Y[6][STATE_SIZE], fp X[6][STATE_SIZE])
{
#pragma HLS INLINE
	fp L[6][6];
	chol6(S,L);

	fp Z[6][STATE_SIZE];
	// forward
	loop_solve1 : for(int i=0;i<6;i++)
		//#pragma HLS PIPELINE II=1
		loop_solve2 : for(int c=0;c<STATE_SIZE;c++)
		{

			fp sum=0;
			loop_solve3 : for(int k=0;k<i;k++)
				sum+=L[i][k]*Z[k][c];
			Z[i][c]=(Y[i][c]-sum)/L[i][i];
		}

	// backward
	loop_solve4 : for(int i=5;i>=0;i--)
		//#pragma HLS PIPELINE II=1
		loop_solve5 : for(int c=0;c<STATE_SIZE;c++)
		{
			fp sum=0;
			loop_solve6 : for(int k=i+1;k<6;k++)
				sum+=L[k][i]*X[k][c];
			X[i][c]=(Z[i][c]-sum)/L[i][i];
		}
}
//===========================================================================//
static void chol3(fp A[3][3], fp L[3][3])
{
	//#pragma HLS INLINE
	loop_chol : for(int i=0;i<3;i++)
	{
		//#pragma HLS PIPELINE II=1
		loop_cho2 : for(int j=0;j<=i;j++)
		{
			fp sum=0;
			loop_cho3 : for(int k=0;k<j;k++)

				sum+=L[i][k]*L[j][k];

			if(i==j)
				L[i][j]=sqrtf(A[i][i]-sum);
			else
				L[i][j]=(A[i][j]-sum)/L[j][j];
		}
		loop_cho4 : for(int j=i+1;j<3;j++)
			L[i][j]=0;
	}
}
//===========================================================================//
static void chol_solve3(fp S[3][3], fp Y[3][STATE_SIZE], fp X[3][STATE_SIZE])
{
#pragma HLS INLINE
	fp L[3][3];
	chol3(S,L);

	fp Z[3][STATE_SIZE];
	// forward
	loop_solve1 : for(int i=0;i<3;i++)
		//#pragma HLS PIPELINE II=1
		loop_solve2 : for(int c=0;c<STATE_SIZE;c++)
		{

			fp sum=0;
			loop_solve3 : for(int k=0;k<i;k++)
				sum+=L[i][k]*Z[k][c];
			Z[i][c]=(Y[i][c]-sum)/L[i][i];
		}

	// backward
	loop_solve4 : for(int i=2;i>=0;i--)
		//#pragma HLS PIPELINE II=1
		loop_solve5 : for(int c=0;c<STATE_SIZE;c++)
		{
			fp sum=0;
			loop_solve6 : for(int k=i+1;k<3;k++)
				sum+=L[k][i]*X[k][c];
			X[i][c]=(Z[i][c]-sum)/L[i][i];
		}
}
//===========================================================================//
void analytic_mag_jacobian
(
		fp state[22],
		fp H_mag[3][22]
)
{
	fp q[4] = {state[0],state[1],state[2],state[3]};
	fp m_ned[3] = {state[16],state[17],state[18]};

	fp Rb2n[3][3];
	quat_to_rot(q,Rb2n);

	fp R_ned2body[3][3];
	lj1 : for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			R_ned2body[i][j] = Rb2n[j][i];

	fp q0 = q[0];
	fp q1 = q[1];
	fp q2 = q[2];
	fp q3 = q[3];
	//
	fp dR_q0[3][3] = {{ 2*q0,-2*q3,2*q2},{2*q3, 2*q0,-2*q1},{-2*q2,2*q1, 2*q0}};
	fp dR_q1[3][3] = {{ 2*q1, 2*q2,2*q3},{2*q2,-2*q1,-2*q0},{ 2*q3,2*q0,-2*q1}};
	fp dR_q2[3][3] = {{-2*q2, 2*q1,2*q0},{2*q1, 2*q2, 2*q3},{-2*q0,2*q3,-2*q2}};
	fp dR_q3[3][3] = {{-2*q3,-2*q0,2*q1},{2*q0,-2*q3, 2*q2},{ 2*q1,2*q2, 2*q3}};
	//
	fp col1[3]= {0};
	fp col2[3]= {0};
	fp col3[3]= {0};
	fp col4[3]= {0};
	lj2 : for(int i=0;i<3;i++)
#pragma HLS PIPELINE II=6
		for(int j=0;j<3;j++)
		{
			col1[i]+=dR_q0[j][i]*m_ned[j];
			col2[i]+=dR_q1[j][i]*m_ned[j];
			col3[i]+=dR_q2[j][i]*m_ned[j];
			col4[i]+=dR_q3[j][i]*m_ned[j];
		}
	//
	fp H[3][22] = {0};
	lj3 : for(int i=0;i<3;i++)
	{
#pragma HLS PIPELINE II=4
		H[i][0] = col1[i];
		H[i][1] = col2[i];
		H[i][2] = col3[i];
		H[i][3] = col4[i];

		H[i][16] = R_ned2body[i][0];
		H[i][17] = R_ned2body[i][1];
		H[i][18] = R_ned2body[i][2];
	}

	H[0][19] = 1;
	H[1][20] = 1;
	H[2][21] = 1;

	lj4 : for(int i=0;i<3;i++)
		for(int j=0;j<22;j++)H_mag[i][j] = H[i][j];
}
//===========================================================================//
void gps_mag_update
(
		fp state[22],
		fp P[22][22],
		fp gps_pos[3],
		fp gps_vel[3],
		fp mag_last[3]
)
{
	fp H_gps[6][STATE_SIZE]={0};
#pragma HLS ARRAY_PARTITION variable=H_gps complete dim = 2

	lg1 : for(int i=0;i<3;i++)
	{
#pragma HLS PIPELINE II=1
		H_gps[i][4+i]=1;
		H_gps[i+3][7+i]=1;
	}

	fp S[6][6]={0};
	lg2 : for(int i=0;i<6;i++)
		for(int j=0;j<6;j++)
			for(int k=0;k<STATE_SIZE;k++)
#pragma HLS PIPELINE II=1
					S[i][j]+=H_gps[i][k]*P[k][k]*H_gps[j][k];

	fp R_gps[6] = {5.169,5.169,5.169,0.0051,0.0051,0.0051};

	lg3 : for(int i=0;i<6;i++)
#pragma HLS PIPELINE II=3
			S[i][i] = S[i][i] + R_gps[i];

	fp Y[6][STATE_SIZE]={0};
#pragma HLS ARRAY_PARTITION variable=Y complete dim = 2
	lg4 : for(int i=0;i<6;i++)
#pragma HLS PIPELINE II=1
		for(int j=0;j<STATE_SIZE;j++)
				Y[i][j]+=H_gps[i][j]*P[j][j];

	fp Kt[6][STATE_SIZE];
#pragma HLS ARRAY_PARTITION variable=Kt complete dim = 2
	chol_solve6(S,Y,Kt);

	fp K[STATE_SIZE][6];
#pragma HLS ARRAY_PARTITION variable=K complete dim = 2
	lg5 : for(int i=0;i<STATE_SIZE;i++)
		for(int j=0;j<6;j++)
			K[i][j]=Kt[j][i];

	fp w3[22][22]={0};
	lg6 : for(int i=0;i<22;i++)
#pragma HLS PIPELINE II=3
		for(int j=0;j<22;j++)
			for(int k=0;k<6;k++)
				w3[i][j]+=K[i][k]*H_gps[k][j];

	lg7 : for(int i=0;i<22;i++)
#pragma HLS PIPELINE II=3
			w3[i][i] = 1 - w3[i][i];

	fp temp[22][22] = {0};
	lg8 : for (int i = 0; i < 22; i++)
	{
	    for (int k = 0; k < 22; k++)
	    {
	        temp[i][k] = w3[i][k] * P[k][k];
	    }
	}

	fp w4[22][22] = {0};
	lg9 : for (int i = 0; i < 22; i++)
	{
//#pragma HLS PIPELINE II=2
		lg91 : for (int j = 0; j < 22; j++)
	    {
//#pragma HLS PIPELINE II=1
			lg92 : for (int k = 0; k < 22; k++)
	        {
#pragma HLS PIPELINE II=1
	            w4[i][j] += temp[i][k] * w3[j][k];
	        }
	    }
	}

	fp w5[22][22]={0};
	lg10 : for(int i=0;i<22;i++)
		for(int j=0;j<22;j++)
			for(int k=0;k<6;k++)
						w5[i][j]+=K[i][k]*R_gps[k]*K[j][k];

	lg11 : for(int i=0;i<22;i++)
		for(int j=0;j<22;j++)
			P[i][j] = w4[i][j] + w5[i][j];

	static fp refloc[3] = {42.2825 ,-72.3430 ,53.0352};

	fp lat0 = refloc[0]*0.0174532925199433;
	fp lon0 = refloc[1]*0.0174532925199433;
	fp alt0 = refloc[2];

	fp lat = gps_pos[0]*0.0174532925199433;
	fp lon = gps_pos[1]*0.0174532925199433;
	fp alt = gps_pos[2];

	fp dy = (lat - lat0) * 6371000;
	fp dx = (lon - lon0) * 4713499;
	fp dz = alt0 - alt;

	fp nedPos_meas[3] = {dy,dx,dz};

	fp z_gps[6];
	z_gps[0] = nedPos_meas[0] - state[4];
	z_gps[1] = nedPos_meas[1] - state[5];
	z_gps[2] = nedPos_meas[2] - state[6];

	z_gps[3] = gps_vel[0] - state[7];
	z_gps[4] = gps_vel[1] - state[8];
	z_gps[5] = gps_vel[2] - state[9];

	fp w6[22] = {0};
	lg12 : for(int i=0;i<22;i++)
		for(int j=0;j<6;j++)
			w6[i]+=K[i][j]*z_gps[j];

	lg13 : for(int i=0;i<22;i++)state[i] = state[i] + w6[i];


	fp H_mag[3][22];
	analytic_mag_jacobian(state,H_mag);

	fp w7[3][22]={0};
	lg14 : for(int i=0;i<3;i++)
		for(int j=0;j<22;j++)
				w7[i][j]+=H_mag[i][j]*P[j][j];

	fp w8[3][3]={0};
	lg15 : for(int i=0;i<3;i++)
//#pragma HLS PIPELINE II=2
		for(int j=0;j<3;j++)
#pragma HLS PIPELINE II=1
			for(int k=0;k<22;k++)
				w8[i][j]+=w7[i][k]*H_mag[j][k];

	fp Rmag = 0.0862;
	fp SS[3][3] = {0};
	lg16 : for(int i=0;i<3;i++)
			SS[i][i] = w8[i][i] + Rmag;

	fp XX[22][3]={0};
	lg17 :for(int i=0;i<22;i++)
#pragma HLS PIPELINE II=3
		for(int j=0;j<3;j++)
				XX[i][j]+=P[i][i]*H_mag[j][i];

	fp YY[3][22]={0};
	lg18 : for(int i=0;i<3;i++)
		for(int j=0;j<22;j++)
			YY[i][j] = XX[j][i];

	fp Ktt[3][22];
	chol_solve3(SS,YY,Ktt);
	fp KK[22][3];
	lg19 : for(int i=0;i<22;i++)
		for(int j=0;j<3;j++)
			KK[i][j] = Ktt[j][i];

	fp q[4] = {state[0],state[1],state[2],state[3]};
	fp m_ned[3] = {state[16],state[17],state[18]};
	fp b_mag[3] = {state[19],state[20],state[21]};

	fp Rb2n[3][3] = {0};
	quat_to_rot(q,Rb2n);
	fp R_ned2body[3][3] = {0};
	lg20 : for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			R_ned2body[i][j] = Rb2n[j][i];

	fp pred_mag[3]={0};
	lg21 : for(int i=0;i<3;i++)
#pragma HLS PIPELINE II=6
		for(int j=0;j<3;j++)
			pred_mag[i]+=R_ned2body[i][j]*m_ned[j];

	for(int i=0;i<3;i++)pred_mag[i] = pred_mag[i] + b_mag[i];

	fp z_mag[3] = {0};
	for(int i=0;i<3;i++)z_mag[i] = mag_last[i] - pred_mag[i];

	fp w9[22]={0};
	for(int i=0;i<22;i++)
		for(int j=0;j<3;j++)
			w9[i]+=KK[i][j]*z_mag[j];

	for(int i=0;i<22;i++)state[i] = state[i] + w9[i];
	//P = (eye(22,'single') - KK * H_mag) * P * (eye(22,'single') - KK * H_mag)' + KK * (Rmag * eye(3,'single')) * KK';
	quat_normalize(&state[0]);
}
//===========================================================================//
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
)
{
	static fp Ts = 0.00625;
	static fp g_ned[3] = {0, 0,9.80665};
	static fp Rpos = 5.169;
	static fp Rvel = 0.0051;
	//===========================================================================//
	IMU_LOOP : for(int k=0;k<IMU_SAMPLES;k++)
		imu_predict(k,state,P,accel_block,gyro_block,orient_block,pos_block);
	gps_mag_update(state,P,gps_pos,gps_vel,mag_last);
}
