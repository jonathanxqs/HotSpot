#ifndef _DEF_H
#include "def.h"
#endif

void intersection_new(double t_1, double r_1, double theta_1, 
					  double phi_1, double t_2, double r_2, 
					  double theta_2, double phi_2, double x_eq[])
{
    double x_1, y_1, z_1;
    double x_2, y_2, z_2;
    double t, x, y, z;
    double ratio;
    double vector[4];

    x_1 = r_1*sin(theta_1)*cos(phi_1);
    y_1 = r_1*sin(theta_1)*sin(phi_1);
    z_1 = r_1*cos(theta_1);
    x_2 = r_2*sin(theta_2)*cos(phi_2);
    y_2 = r_2*sin(theta_2)*sin(phi_2);
    z_2 = r_2*cos(theta_2);

    vector[0] = t_1 - t_2;
    vector[1] = x_1 - x_2;
    vector[2] = y_1 - y_2;
    vector[3] = z_1 - z_2;
    ratio = abs(z_2)/vector[3];

    t = t_2 + ratio*vector[0];
    x = x_2 + ratio*vector[1];
    y = y_2 + ratio*vector[2];
    z = z_2 + ratio*vector[3];

    x_eq[0] = t;
    x_eq[1] = sqrt(x*x + y*y);
    x_eq[2] = Pi/2;
    x_eq[3] = asin(y/x_eq[1]);

    if (x < 0)  x_eq[3] = Pi - x_eq[3];

    return;
}
