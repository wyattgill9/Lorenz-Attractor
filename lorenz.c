#include<stdio.h>
#include<time.h>

#define SIGMA 10.0
#define RHO 28.0
#define BETA (8.0/3.0)

void lorenz_rk4(double *x, double *y, double *z, double h, int steps, FILE *datafile) {
    double x_val = *x;
    double y_val = *y;
    double z_val = *z;
    
    for (int i = 0; i < steps; i++) {
        double dx1 = SIGMA * (y_val - x_val);
        double dy1 = x_val * (RHO - z_val) - y_val;
        double dz1 = x_val * y_val - BETA * z_val;
        
        double x1 = x_val + h * dx1/2.0;
        double y1 = y_val + h * dy1/2.0;
        double z1 = z_val + h * dz1/2.0;
        
        double dx2 = SIGMA * (y1 - x1);
        double dy2 = x1 * (RHO - z1) - y1;
        double dz2 = x1 * y1 - BETA * z1;
        
        double x2 = x_val + h * dx2/2.0; 
        double y2 = y_val + h * dy2/2.0; 
        double z2 = z_val + h * dz2/2.0; 
        
        double dx3 = SIGMA * (y2 - x2);
        double dy3 = x2 * (RHO - z2) - y2;
        double dz3 = x2 * y2 - BETA * z2;
        
        double x3 = x_val + h * dx3;
        double y3 = y_val + h * dy3;
        double z3 = z_val + h * dz3;
        
        double dx4 = SIGMA * (y3 - x3);
        double dy4 = x3 * (RHO - z3) - y3;
        double dz4 = x3 * y3 - BETA * z3;
        
        x_val += h * (dx1 + 2*dx2 + 2*dx3 + dx4)/6.0;
        y_val += h * (dy1 + 2*dy2 + 2*dy3 + dy4)/6.0;
        z_val += h * (dz1 + 2*dz2 + 2*dz3 + dz4)/6.0;
        
        fprintf(datafile, "%d %lf %lf %lf\n", i, x_val, y_val, z_val);
    }
    
    *x = x_val;
    *y = y_val;
    *z = z_val;
}

int main() {
    double x = 1.0, y = 1.0, z = 1.0;
    double h = 0.01;
    int steps = 10000;
    
    FILE *datafile = fopen("lorenz_data.dat", "w");
    if (datafile == NULL) {
        printf("Error opening file!\n");
        return 1;
    }
    
    lorenz_rk4(&x, &y, &z, h, steps, datafile);
    fclose(datafile);
    
    FILE *gnuplotScript = fopen("plot_lorenz.gp", "w");
    if (gnuplotScript == NULL) {
        printf("Error creating gnuplot script!\n");
        return 1;
    }
    
    fprintf(gnuplotScript, "set terminal png size 1200,900\n");
    fprintf(gnuplotScript, "set output 'lorenz_attractor.png'\n");
    fprintf(gnuplotScript, "set title 'Lorenz Attractor'\n");
    fprintf(gnuplotScript, "set xlabel 'X'\n");
    fprintf(gnuplotScript, "set ylabel 'Y'\n");
    fprintf(gnuplotScript, "set zlabel 'Z'\n");
    fprintf(gnuplotScript, "splot 'lorenz_data.dat' using 2:3:4 with lines lc rgb 'blue' title 'Trajectory'\n");
    fclose(gnuplotScript);
    
    printf("Data saved to lorenz_data.dat\n");
    printf("To generate the plot, run: gnuplot plot_lorenz.gp\n");
    
    return 0;
}
