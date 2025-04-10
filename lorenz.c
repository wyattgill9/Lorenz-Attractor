#include<stdio.h>
#include<time.h>

#define SIGMA 10.0
#define RHO 28.0
#define BETA (8.0/3.0)

void lorenz(double x, double y, double z, double h, int steps) {
    for (int i = 0; i < steps; i++) {
        double dx = SIGMA * (y - x);
        double dy = x * (RHO - z) - y;
        double dz = x * y - BETA * z;

        x += h * dx;
        y += h * dy;
        z += h * dz;
        
        // printf("%lf %lf %lf\n", x, y, z);
    }
}

void lorenz_rk4(double x, double y, double z, double h, int steps) {
    for (int i = 0; i < steps; i++) {
        double dx1 = SIGMA * (y - x);
        double dy1 = x * (RHO - z) - y;
        double dz1 = x * y - BETA * z;

        double x1 = x + h * dx1/2.0;
        double y1 = y + h * dy1/2.0;
        double z1 = z + h * dz1/2.0;

        double dx2 = SIGMA * (y1 - x1);
        double dy2 = x1 * (RHO - z1) - y1;
        double dz2 = x1 * y1 - BETA * z1;

        double x2 = x + h * dx2/2.0; 
        double y2 = y + h * dy2/2.0; 
        double z2 = z + h * dz2/2.0; 

        double dx3 = SIGMA * (y2 - x2);
        double dy3 = x2 * (RHO - z2) - y2;
        double dz3 = x2 * y2 - BETA * z2;

        double x3 = x + h * dx3;
        double y3 = y + h * dy3;
        double z3 = z + h * dz3;

        double dx4 = SIGMA * (y3 - x3);
        double dy4 = x3 * (RHO - z3) - y3;
        double dz4 = x3 * y3 - BETA * z3;

        x += h * (dx1 + 2*dx2 + 2*dx3 + dx4)/6.0;
        y += h * (dy1 + 2*dy2 + 2*dy3 + dy4)/6.0;
        z += h * (dz1 + 2*dz2 + 2*dz3 + dz4)/6.0;

        // printf("%lf %lf %lf\n", x, y, z);
    }
}

int main() {
    clock_t start_time, end_time;

    double x = 1.0, y = 1.0, z = 1.0;
    double h = 0.01;
    int steps = 100;

    start_time = clock();
    // lorenz(x, y, z, h, steps);
    lorenz_rk4(x, y, z, h, steps);
    end_time = clock();

    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Finished in %.7f seconds\n", elapsed_time);

    /*
    BENCHMARKS
        rk4 -        Finished in 0.0000120 seconds
        rk4 + O3 -   Finished in 0.0000050 seconds
        rk4 + O3 + native - Finished in 0.0000030 seconds
    
        naive -      Finished in 0.0000070 seconds
        naive + 03 - Finished in 0.0000050 seconds

    Conclusion:
        Equal time for both w/t compiler optimizations so its better to use rk4 anyways
    
    */

    return 0;
}
