#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


void leer_red(const char* filename, int** vecinos, int** grados, int* total_nodos) {



    // Esta función recibe como entrada el nombre de un archivo de texto donde se encuentra representada la red. Te saca dos vecotres unidimensionales:
    // La componente i-esima de grados nos dice el grado del nodo i-esimo de la red
    // vecinos guarda secuencialmente los vecinos del archivo de la matriz de adyacencia

    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error al abrir el archivo");
        exit(EXIT_FAILURE);
    }

    int max_node = -1;
    int a, b;
    int edge_count = 0;

    // Primera pasada: encontrar el máximo nodo y contar aristas
    while (fscanf(file, "%d %d", &a, &b) == 2) {
        if (a > max_node) max_node = a;
        if (b > max_node) max_node = b;
        edge_count++;
    }

    *total_nodos = max_node + 1;
    rewind(file);

    // Reservar memoria para grados y vecinos
    *grados = calloc(*total_nodos, sizeof(int));
    *vecinos = malloc(edge_count * 2 * sizeof(int));

    // Segunda pasada: contar grados
    int i=0;
    while (fscanf(file, "%d %d", &a, &b) == 2) {
        (*grados)[a]++;
        (*vecinos)[i]=b;
        i++;
    }
    fclose(file);
}


double sumando_i(double K, double betta, double delta, int* vecinos, int grado, double* x, int i) {
    double* w_j = malloc(grado * sizeof(double));  // array dinámico para pesos
    if (w_j == NULL) {
        fprintf(stderr, "Error: no se pudo asignar memoria\n");
        exit(1);
    }

    double denominador = 0.0;
    int l;
    double sumando = 0.0;

    // Calcular denominador
    for (l = 0; l < grado; l++) {
        denominador +=pow( fabs((x[i] - x[vecinos[l]])) + delta,-betta);
    }

    // Calcular pesos w_j
    for (l = 0; l < grado; l++) {
        double numerador = pow(fabs((x[i] - x[vecinos[l]])) + delta,-betta);
        w_j[l] = numerador/denominador;
    }

    // Calcular sumando
    for (l = 0; l < grado; l++) {
        sumando += w_j[l] * tanh((double)x[vecinos[l]]);
    }

    free(w_j);  // liberar memoria asignada
    return K * sumando;
}

void derivada(double K, double betta, double delta, int* vecinos, int *grados, double* x, double* derivada_x, int total_nodos){
int indice=0;
for(int i=0;i<total_nodos;i++){
    derivada_x[i]=-x[i]+sumando_i(K,betta,delta,vecinos+indice,grados[i],x,i);
    indice=indice+grados[i];
}

}

void rk4_step(double* x, int total_nodos, double dt,
              double K, double betta, double delta,
              int* vecinos, int* grados) {

    double* k1 = malloc(total_nodos * sizeof(double));
    double* k2 = malloc(total_nodos * sizeof(double));
    double* k3 = malloc(total_nodos * sizeof(double));
    double* k4 = malloc(total_nodos * sizeof(double));
    double* x_temp = malloc(total_nodos * sizeof(double));

    if (!k1 || !k2 || !k3 || !k4 || !x_temp) {
        fprintf(stderr, "Error: no se pudo asignar memoria\n");
        exit(1);
    }

    // k1 = f(x)
    derivada(K, betta, delta, vecinos, grados, x, k1, total_nodos);

    // k2 = f(x + dt/2 * k1)
    for (int i = 0; i < total_nodos; i++)
        x_temp[i] = x[i] + 0.5 * dt * k1[i];
    derivada(K, betta, delta, vecinos, grados, x_temp, k2, total_nodos);

    // k3 = f(x + dt/2 * k2)
    for (int i = 0; i < total_nodos; i++)
        x_temp[i] = x[i] + 0.5 * dt * k2[i];
    derivada(K, betta, delta, vecinos, grados, x_temp, k3, total_nodos);

    // k4 = f(x + dt * k3)
    for (int i = 0; i < total_nodos; i++)
        x_temp[i] = x[i] + dt * k3[i];
    derivada(K, betta, delta, vecinos, grados, x_temp, k4, total_nodos);

    // Actualizar x usando combinación ponderada
    for (int i = 0; i < total_nodos; i++) {
        x[i] += dt * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) / 6.0;
    }

    // Liberar memoria
    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(x_temp);
}

double condiciones_iniciales(double K, int total_nodos, double *x0){
double gamma;
if(K<1){
    for(int j=0;j<total_nodos;j++){
        x0[j]=-1 + ((double)rand() / RAND_MAX)*2;
    }
            gamma=0.002;
}
    else{

        for(int j=0;j<total_nodos;j++){
        x0[j]=-K + ((double)rand() / RAND_MAX)*2*K;
    }
    gamma=0.002*K;
}
return gamma;
}


void polarizacion(double *x, int total_nodos, double *abs_op_media, double *desvest){
double promedio=0;
    for(int i=0;i<total_nodos;i++){
        promedio=promedio+x[i];
    }
    promedio=promedio/total_nodos;
    *abs_op_media= fabs(promedio);
    *desvest=0;
        for(int i=0;i<total_nodos;i++){
        *desvest=*desvest+(x[i]-promedio)*(x[i]-promedio);
    }
    *desvest=sqrt(*desvest/total_nodos);
}

void copia_conf(double* x_original, double *x_copiado, int total_nodos){
for(int i=0;i<total_nodos;i++){
    x_copiado[i]=x_original[i];
}

}
void escribe_evolucion_polarizacion(double tiempo, double x_med, double x_desvest, const char* filename_output){
    FILE *archivo = fopen(filename_output, "a");  // "a" para añadir sin sobrescribir

    if (archivo == NULL) {
        fprintf(stderr, "Error al abrir el archivo %s\n", filename_output);
        return;
    }

    fprintf(archivo, "%lf\t%lf\t%lf\n", tiempo, x_med, x_desvest);
    fclose(archivo);
}



void simula_polarizacion(const char* filename_input, int N_pasos, double dt,double K, double betta, const char* filename_output){
    FILE *archivo = fopen(filename_output, "w"); //Para crear el archivo
    fprintf(archivo, "%lf\t%lf\t%d\t%lf\n", K, betta, N_pasos, dt);

    fclose(archivo);
    int* vecinos;
    int* grados;
    int total_nodos = 0;
    double delta;
    double x_med, x_desvest;
    leer_red(filename_input, &vecinos, &grados, &total_nodos);
    double x[total_nodos];
    delta= condiciones_iniciales(K,total_nodos,x);
    polarizacion(x,total_nodos, &x_med, &x_desvest);
    escribe_evolucion_polarizacion(0,x_med, x_desvest,filename_output);
    for(int j=0;j<N_pasos;j++){
        rk4_step(x,total_nodos,dt,K, betta,delta,vecinos, grados);
        polarizacion(x,total_nodos, &x_med, &x_desvest);
        escribe_evolucion_polarizacion(dt*(j+1),x_med, x_desvest,filename_output);
    }

}



int main() {
    srand(time(NULL)); // Inicializa aleatoriedad

    const char* archivo_red = "C:/Users/USUARIO/Downloads/CAOS/CaOtIcOs/ARCHIVOS_REDES/ER/ER_0.txt";
    const char* archivo_salida = "salida_polarizacion.txt";

    double dt = 0.1;
    int N_pasos = 100;
    double K = 0.5;
    double betta = 1;
    int *vecinos;
    int *grados;
    int total_nodos=0;
    simula_polarizacion(archivo_red, N_pasos, dt, K, betta, archivo_salida);

    return 0;
}
