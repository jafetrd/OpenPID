/* pid.h - PID Controller Header */
#ifndef PID_H
#define PID_H

typedef struct {
    float kp;  // Proportional gain
    float ki;  // Integral gain
    float kd;  // Derivative gain
    float prev_error;  // Previous error
    float integral;  // Integral term
} PIDController;

void pid_init(PIDController *pid, float kp, float ki, float kd);
float pid_compute(PIDController *pid, float setpoint, float measurement, float dt);

#endif /* PID_H */

/* pid.c - PID Controller Implementation */
#include "pid.h"

void pid_init(PIDController *pid, float kp, float ki, float kd) {
    pid->kp = kp;
    pid->ki = ki;
    pid->kd = kd;
    pid->prev_error = 0.0f;
    pid->integral = 0.0f;
}

float pid_compute(PIDController *pid, float setpoint, float measurement, float dt) {
    float error = setpoint - measurement;
    pid->integral += error * dt;
    float derivative = (error - pid->prev_error) / dt;
    float output = (pid->kp * error) + (pid->ki * pid->integral) + (pid->kd * derivative);
    pid->prev_error = error;
    return output;
}
