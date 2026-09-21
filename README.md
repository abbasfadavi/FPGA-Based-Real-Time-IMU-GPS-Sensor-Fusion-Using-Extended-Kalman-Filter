# FPGA-Based Real-Time IMU/GPS Sensor Fusion Using Extended Kalman Filter

A hardware-oriented implementation of an IMU/GPS sensor fusion system based on the Extended Kalman Filter (EKF), developed for FPGA deployment using High-Level Synthesis (HLS).

The main goal of this project is to demonstrate how a MATLAB-based navigation algorithm can be transformed into a synthesizable FPGA architecture while preserving numerical behavior and achieving deterministic real-time execution.

The design was implemented using Vitis HLS 2023.1 and synthesized for a Xilinx Kintex-7 XC7K410T-FFG900-2 FPGA.

---

## Overview

Inertial Measurement Units (IMUs) provide high-rate acceleration and angular-rate measurements but suffer from accumulated drift. GPS provides absolute position and velocity measurements but typically operates at a significantly lower update rate.

This project combines both sensors using an Extended Kalman Filter.

Target configuration:
- IMU update rate: 160 Hz
- GPS update rate: 1 Hz
- 160 IMU prediction steps per GPS correction cycle

---

## Motivation

A MATLAB implementation of an EKF cannot generally be transferred directly to FPGA hardware. Several software-oriented constructs are unsuitable or inefficient for HLS, including dynamically sized arrays, toolbox-specific objects, data-dependent loops, iterative Jacobian computation, general matrix inversion, dynamically controlled execution, and software-oriented quaternion operations.

This project restructures the original MATLAB implementation into a hardware-compatible form.

---

## State Vector

The original MATLAB navigation model used a 22-state EKF. The FPGA-oriented implementation reduces the state vector to 16 states:

Position                 3
Velocity                 3
Quaternion               4
Accelerometer bias       3
Gyroscope bias           3
--------------------------------
Total                   16

The geomagnetic field and magnetometer bias states were removed. Consequently, the covariance matrix was reduced from 22×22 to 16×16, reducing memory requirements and computational complexity.

---

## EKF Processing Flow

The implementation consists of two main processing stages.

### IMU Prediction

The prediction stage processes accelerometer and gyroscope measurements and updates:
- position
- velocity
- orientation quaternion
- sensor biases
- covariance matrix

### GPS Update

When a GPS measurement becomes available, the EKF performs a correction using GPS position and velocity.

Measurement vector : 6×1
Observation matrix : 6×16
Kalman gain        : 16×6

---

## MATLAB-to-HLS Conversion

Several modifications were required to convert the MATLAB navigation implementation into synthesizable HLS C/C++.

### 1. Removal of Toolbox Dependencies

MATLAB-specific functionality was replaced with synthesizable implementations. For example, lla2ned() was replaced by a lightweight coordinate transformation based on a spherical Earth approximation. Reference-dependent constants are calculated once during initialization instead of repeatedly during runtime.

### 2. Quaternion Operations

MATLAB quaternion objects were replaced by fixed-size arrays. Quaternion multiplication was implemented explicitly using the Hamilton product.

### 3. Static Loop Bounds

The original MATLAB implementation used control structures that were not appropriate for hardware synthesis. They were replaced by statically bounded loops:

50 GPS blocks
    └── 160 IMU samples per GPS block

Static bounds enable pipelining, loop analysis, scheduling, and resource allocation.

---

## Hardware-Oriented Optimizations

### Analytical Jacobian

Iterative finite-difference Jacobian evaluation was replaced by an analytical formulation, reducing repeated nonlinear evaluations and computational overhead.

### Cholesky-Based Kalman Gain

Direct matrix inversion is expensive in FPGA hardware. Kalman gain computation was therefore reformulated using Cholesky factorization for improved numerical stability and more hardware-friendly implementation.

### Covariance Update

A simplified covariance update was used instead of the full Joseph-form covariance update, reducing the number of matrix multiplications required by the EKF correction stage.

---

## HLS Optimizations

Several HLS directives were used to improve parallelism and memory bandwidth, including:

#pragma HLS PIPELINE
#pragma HLS ARRAY_PARTITION

Main optimization targets:
- IMU prediction loop
- matrix multiplication
- covariance propagation
- Kalman gain computation

---

## Data Structures

State vector          : 16×1   — Position, velocity, quaternion and IMU biases
Covariance matrix     : 16×16  — State estimation uncertainty
Process noise matrix  : 16×16  — IMU process noise
GPS noise matrix      : 6×6    — GPS measurement noise
Observation matrix    : 6×16   — GPS measurement model
Kalman gain           : 16×6   — EKF correction gain

All arrays use compile-time fixed dimensions to maintain HLS synthesizability.

---

## FPGA Target

Tool        : Vitis HLS 2023.1
FPGA        : Xilinx Kintex-7
Device      : XC7K410T-FFG900-2
Target Clock: 20 ns
Target Freq.: 50 MHz

HLS-estimated clock period: 16.838 ns
Estimated maximum operating frequency: approximately 59.4 MHz

---

## Performance

Complete EKF update latency:
55,922–60,266 clock cycles

At 50 MHz:
1.118–1.205 ms per estimation cycle

Individual stages:
- IMU prediction: 187 cycles (~3.74 µs)
- GPS update: 25,678–30,022 cycles (~0.514–0.600 ms)

The GPS correction stage dominates execution time because of the matrix operations required for covariance and Kalman gain calculations.

---

## FPGA Resource Utilization

LUT       : 44,317 (17%)
Flip-Flop : 42,955 (8%)
DSP       : 191 (12%)
BRAM      : 37 (2%)

The majority of resources are consumed by the GPS update stage and its floating-point matrix operations.

---

## Numerical Verification

The FPGA/HLS implementation was validated against the MATLAB reference implementation using identical input sequences.

Compared outputs include:
- estimated position
- estimated velocity
- orientation quaternion
- covariance matrix

The observed numerical difference remained below approximately 10^-3% for the evaluated dataset, demonstrating close numerical agreement between the MATLAB reference model and the HLS implementation.

---

## Floating-Point Implementation

The current implementation retains double-precision floating-point operations for several nonlinear functions, including:
- sqrt
- sin
- cos
- norm

This choice was intentionally made to simplify functional verification against the MATLAB reference implementation.

Therefore, this repository should primarily be considered a functional and architectural FPGA/HLS implementation rather than a fully resource-minimized production implementation.

A future fixed-point version could significantly reduce DSP usage, LUT usage, latency, and power consumption. Possible replacements include fixed-point arithmetic, CORDIC, lookup tables, and mixed-precision arithmetic.

---

## Dataset

Functional verification was performed using IMU/GPS data derived from the MATLAB Navigation Toolbox IMU/GPS navigation example. The same input sequence was applied to both the MATLAB reference implementation and the HLS implementation to enable direct numerical comparison.

---

## Main Engineering Challenges

1. MATLAB toolbox dependency removal
2. Dynamically sized data structure elimination
3. Fixed-size matrix design
4. Synthesizable quaternion arithmetic
5. Static control-flow generation
6. Analytical Jacobian implementation
7. Hardware-efficient matrix computation
8. HLS pipelining
9. Memory partitioning
10. MATLAB/HLS numerical validation

---

## Possible Future Work

- complete fixed-point conversion
- automatic word-length optimization
- mixed-precision EKF implementation
- CORDIC-based nonlinear functions
- deeper matrix-operation pipelining
- parallel matrix multiplication
- resource/latency design-space exploration
- comparison between HLS and handwritten HDL
- power consumption evaluation
- deployment on the physical Kintex-7 platform
- integration with a complete real-time navigation system

---

## Intended Use

This repository is intended primarily for:
- FPGA engineers
- DSP engineers
- navigation-system developers
- HLS researchers
- students working on IMU/GPS fusion
- engineers interested in MATLAB-to-FPGA migration

The project may also serve as a reference implementation for studying the transformation of matrix-intensive estimation algorithms from software to FPGA hardware.

---

## Tools

- MATLAB
- MATLAB Navigation Toolbox
- Vitis HLS 2023.1
- Xilinx Kintex-7 FPGA

---

## Author

Abbas Fadavi
FPGA / Real-Time DSP Engineer

GitHub:
https://github.com/abbasfadavi

LinkedIn:
https://www.linkedin.com/in/abbasfadavi/

---

## Disclaimer

This project is a research and engineering prototype intended for educational and experimental use.

The current implementation prioritizes functional equivalence with the MATLAB reference model and deterministic FPGA execution. Additional verification, fixed-point optimization, hardware-in-the-loop testing, and application-specific validation would be required before deployment in safety-critical navigation systems.
