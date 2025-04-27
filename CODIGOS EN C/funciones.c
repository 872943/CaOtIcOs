#include "funciones.h"
#include <stdbool.h>
#include <dirent.h>


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
// Calcula el primer sumando de la EDO 1 del articulo

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
// Esta funcion calcula la derivada de la EDO 1 del articulo
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

double condiciones_iniciales(double K, int total_nodos, double *x0) {
    if (total_nodos <= 0) {
        fprintf(stderr, "Error: total_nodos debe ser mayor que 0\n");
        exit(EXIT_FAILURE);
    }

    double gamma;
    if (K < 1) {
        for (int j = 0; j < total_nodos; j++) {
            // Generar valores entre -1 y 1
            x0[j] = -1.0 + ((double)rand() / (double)RAND_MAX) * 2.0;
        }
        gamma = 0.002;
    } else {
        for (int j = 0; j < total_nodos; j++) {
            // Generar valores entre -K y K
            x0[j] = -K + ((double)rand() / (double)RAND_MAX) * (2.0 * K);
        }
        gamma = 0.002 * K;
    }

    return gamma;
}

void polarizacion(double *x, int total_nodos, double *op_media, double *desvest){
double promedio=0;
    for(int i=0;i<total_nodos;i++){
        promedio=promedio+x[i];
    }
    promedio=promedio/total_nodos;
    *op_media= promedio;
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

double cacula_velocidad_modulo(double K, double betta, double delta, int* vecinos, int *grados, double* x, int total_nodos){ 

    double * velocidad = malloc(total_nodos * sizeof(double));
    if (velocidad == NULL) {
        fprintf(stderr, "Error: no se pudo asignar memoria\n");
        exit(1);
    }
    derivada(K,betta,delta,vecinos,grados,x,velocidad,total_nodos);

    double suma = 0.0;
    for (int i = 0; i < total_nodos; i++) {
        suma += velocidad[i] * velocidad[i];
    }
    return sqrt(suma);
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


// Busca el índice más alto entre los archivos ER_*.txt en la carpeta de salida
int obtener_siguiente_indice(const char* carpeta) {
    DIR* dir;
    struct dirent* entry;
    int max_index = -1;

    dir = opendir(carpeta);
    if (dir == NULL) {
        perror("No se pudo abrir el directorio de salida");
        return 0;
    }

    while ((entry = readdir(dir)) != NULL) {
        int idx;
        if (sscanf(entry->d_name, "ER_%d.txt", &idx) == 1) {
            if (idx > max_index) {
                max_index = idx;
            }
        }
    }

    closedir(dir);
    return max_index + 1;
}

void escribe_evolucion_individual(char*filename_output,int N_nodos, double *x, double tiempo){
    FILE *archivo = fopen(filename_output, "a");  // "a" para añadir sin sobrescribir

    if (archivo == NULL) {
        fprintf(stderr, "Error al abrir el archivo %s\n", filename_output);
        return;
    }

    fprintf(archivo, "%lf\t", tiempo);
    for(int i=0;i<N_nodos;i++){
        fprintf(archivo, "%lf\t", x[i]);
    }
    fprintf(archivo, "\n");
    fclose(archivo);

}

void escribe_velocidad_modulo(char*filename_output,double tiempo, double velocidad){
    FILE *archivo = fopen(filename_output, "a");  // "a" para añadir sin sobrescribir

    if (archivo == NULL) {
        fprintf(stderr, "Error al abrir el archivo %s\n", filename_output);
        return;
    }

    fprintf(archivo, "%lf\t%lf\n", tiempo,velocidad);
    fclose(archivo);

}

void evolucion_persona_a_persona(char*filename_input,char*filename_output, int N_pasos, double dt, double K, double betta, char*filename_output_velocidad){
    // Esta función simula la evolución de la polarización de una red a lo largo del tiempo y guarda los resultados en un archivo de salida.
    // La red se lee desde un archivo de entrada y se simula durante N_pasos pasos de tiempo con un tamaño de paso dt.
    // Los parámetros K y betta son constantes que afectan la dinámica de la red. number_name es un número que se utiliza para identificar la simulación.
    FILE *archivo_velocidad = fopen(filename_output_velocidad, "w"); //Para crear el archivo
    fprintf(archivo_velocidad, "%lf\t%lf\t%d\t%lf\t%s\n", K, betta, N_pasos, dt, filename_input);
    FILE *archivo = fopen(filename_output, "w"); //Para crear el archivo
    fprintf(archivo, "%lf\t%lf\t%d\t%lf\t%s\n", K, betta, N_pasos, dt, filename_input);
    int* vecinos;
    int* grados;
    int total_nodos = 0;
    double desvest;
    double op_media;
    double delta;
    leer_red(filename_input, &vecinos, &grados, &total_nodos);
    double x[total_nodos];
    delta= condiciones_iniciales(K,total_nodos,x);
    escribe_evolucion_individual(filename_output,total_nodos,x,dt*0.0);
    for(int j=0;j<N_pasos;j++){
        rk4_step(x,total_nodos,dt,K, betta,delta,vecinos, grados);
        escribe_evolucion_individual(filename_output,total_nodos,x,dt*(j+1));
        double velocidad_modulo=cacula_velocidad_modulo(K,betta,delta,vecinos,grados,x,total_nodos);
        escribe_velocidad_modulo(filename_output_velocidad,dt*(j+1),velocidad_modulo);
    }

    fclose(archivo_velocidad);
    fclose(archivo);

}

bool esta_termalizada(double K, double betta, double delta, int* vecinos, int* grados, double* x, int total_nodos) {
    double a = 0.03;  // Parámetro para determinar si el sistema está termalizado
    double v = cacula_velocidad_modulo(K, betta, delta, vecinos, grados, x, total_nodos);
    if (K<1){
        return (v<sqrt(total_nodos)*a);
    }else{
        return (v<sqrt(total_nodos)*a*K);

    }
}

void muchas_simulaciones_ER(int N_sim, int N_pasos, double dt, double K, double betta, int number_name) {
    char* direccion_input = "C:\\Users\\HP\\Desktop\\FISICA\\3 (2024-2025)\\segundo cuatri\\caos\\trabajo\\CaOtIcOs\\ARCHIVOS_REDES\\ER";
    char* carpeta_output = "C:\\Users\\HP\\Desktop\\FISICA\\3 (2024-2025)\\segundo cuatri\\caos\\trabajo\\CaOtIcOs\\Resultados (PARTE 4)\\Evolucion temporal (promedio y desvest)\\ER";

    int indice_salida = obtener_siguiente_indice(carpeta_output);
    int indice_salida_inicial = indice_salida;

    for (int j = 0; j < N_sim; j++) {
        char filename_input[512];
        char filename_output[512];

        int id_input = number_name - j;
        int id_output = indice_salida++;
        printf("1");

        sprintf(filename_input, "%s_%d.txt", direccion_input, id_input);
        printf("me cago en %s",filename_input);
        sprintf(filename_output, "%s\\ER_%d.txt", carpeta_output, id_output);

        simula_polarizacion(filename_input, N_pasos, dt, K, betta, filename_output);
    }

    // Guardar historial

    char historial_path[512];

    sprintf(historial_path, "%s\\ER_Historial.txt", carpeta_output);

    FILE* historial = fopen(historial_path, "a");
    if (historial == NULL) {
        perror("No se pudo abrir el archivo de historial");
        return;
    }

    time_t t = time(NULL);
    struct tm tm = *localtime(&t);

    fprintf(historial, "--------------------------------------------------------------------------------\n");
    fprintf(historial, "Fecha: %04d-%02d-%02d %02d:%02d:%02d\n", 
            tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday,
            tm.tm_hour, tm.tm_min, tm.tm_sec);
    fprintf(historial, "Simulaciones: %d\n", N_sim);
    fprintf(historial, "Redes (input): ER_%d.txt a ER_%d.txt\n", number_name, number_name - (N_sim - 1));
    fprintf(historial, "Resultados (output): ER_%d.txt a ER_%d.txt\n", 
            indice_salida_inicial, indice_salida - 1);
    fprintf(historial, "Parámetros -> K: %.2lf | β: %.2lf | dt: %.4lf | N_pasos: %d\n", 
            K, betta, dt, N_pasos);
    fprintf(historial, "--------------------------------------------------------------------------------\n\n");

    fclose(historial);
}

void muchas_simulaciones_WS(int N_sim, int N_pasos, double dt, double K, double betta, int number_name) {
    char* direccion_input = "C:\\Users\\USUARIO\\Downloads\\CAOS\\CaOtIcOs\\ARCHIVOS_REDES\\WS\\WS";
    char* carpeta_output = "C:\\Users\\USUARIO\\Downloads\\CAOS\\CaOtIcOs\\Resultados (Parte 0)\\Evolucion temporal (promedio y desvest)\\WS";

    int indice_salida = obtener_siguiente_indice(carpeta_output);
    int indice_salida_inicial = indice_salida;

    for (int j = 0; j < N_sim; j++) {
        char filename_input[512];
        char filename_output[512];

        int id_input = number_name - j;
        int id_output = indice_salida++;

        sprintf(filename_input, "%s_%d.txt", direccion_input, id_input);
        sprintf(filename_output, "%s\\WS_%d.txt", carpeta_output, id_output);

        simula_polarizacion(filename_input, N_pasos, dt, K, betta, filename_output);
    }

    // Guardar historial
    char historial_path[512];
    sprintf(historial_path, "%s\\WS_Historial.txt", carpeta_output);

    FILE* historial = fopen(historial_path, "a");
    if (historial == NULL) {
        perror("No se pudo abrir el archivo de historial");
        return;
    }

    time_t t = time(NULL);
    struct tm tm = *localtime(&t);

    fprintf(historial, "--------------------------------------------------------------------------------\n");
    fprintf(historial, "Fecha: %04d-%02d-%02d %02d:%02d:%02d\n", 
            tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday,
            tm.tm_hour, tm.tm_min, tm.tm_sec);
    fprintf(historial, "Simulaciones: %d\n", N_sim);
    fprintf(historial, "Redes (input): WS_%d.txt a WS_%d.txt\n", number_name, number_name - (N_sim - 1));
    fprintf(historial, "Resultados (output): WS_%d.txt a WS_%d.txt\n", 
            indice_salida_inicial, indice_salida - 1);
    fprintf(historial, "Parámetros -> K: %.2lf | β: %.2lf | dt: %.4lf | N_pasos: %d\n", 
            K, betta, dt, N_pasos);
    fprintf(historial, "--------------------------------------------------------------------------------\n\n");

    fclose(historial);
}

double desviacion_estandar(double arr[], int tam, double media) {
    // Función que calcula la desviación estándar de un array de doubles
    double suma = 0.0;
    for (int i = 0; i < tam; i++) {
        suma += (arr[i] - media) * (arr[i] - media);
    }
    return sqrt(suma / tam);
}
bool es_polarizada(const char* filename) {
    //le pasas el nombre de un archivo con la evolución temporal de una red y a partir de los
    //últimos valores registrados para la media y la varianza, te dice si es polarizada (true) o no (false)
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Error al abrir el archivo %s\n", filename);
        exit(1);
    }

    double tiempo, x_medio, desviacion;

    double t, x, desv;
    while (fscanf(file, "%lf %lf %lf", &t, &x, &desv) == 3) {
        tiempo = t;
        x_medio = x;
        desviacion = desv;
    }

    fclose(file);

    
    //el criterio que vamos a seguir para decidir si es polarizada o no es el siguiente:
           //si la media=0 y la varianza!=0, es polarizada
           //si la media!=0 y la varianza!=0, entonces sera polarizada si la media es menor
           //para el resto de casos NO es polarizada

    double eps=0.0001; //tolerancia para evitar problemas de redondeo, REVISAR
    
    if (fabs(x) < eps && desv > eps) {
        return true;  // Polarizada
    } else if (fabs(x) > eps && desv > eps) {
        return fabs(x) < desv;  // Polarizada si la media es menor que la desviación
    } else {
        return false;  // No polarizada
    }

}



double valor_medio(double arr[], int tam) {
// Función que calcula el valor medio de un array de doubles
    double suma = 0.0;
    for (int i = 0; i < tam; i++) {
        suma += arr[i];
    }
    return suma / tam;
}

void frac_polarizado(int N_redes, int rede_ini, double K, double betta, double N_pasos, double dt) {



    char carpeta_K[64];
    char carpeta_betta[64];
    char carpeta_output[128];

    // Crear nombres de carpetas a partir de K y betta
    sprintf(carpeta_K, "%.2lf", K);
    sprintf(carpeta_betta, "%.2lf", betta);
    sprintf(carpeta_output, "%s\\%s", carpeta_K, carpeta_betta);

    // Intentar crear carpeta K
    if (mkdir_p(carpeta_K) != 0 && errno != EEXIST) {
        fprintf(stderr, "❌ Error: No se pudo crear la carpeta '%s'.\n", carpeta_K);
        return;
    }

    // Intentar crear carpeta betta
    if (mkdir_p(carpeta_output) != 0 && errno != EEXIST) {
        fprintf(stderr, "❌ Error: No se pudo crear la carpeta '%s'.\n", carpeta_output);
        return;
    }

    int polarizadas = 0, ult;
    int no_polarizadas = 0;
    char filename_input[512];
    double op_media, desvest;
    char filename_output[512];

    ult=obtener_siguiente_indice(carpeta_output);

    for (int i = 0; i < N_redes; i++) {
        printf("+1");
        sprintf(filename_input, "ARCHIVOS_REDES\\ER\\ER_%d.txt", rede_ini + i);
        sprintf(filename_output, "%s\\ER_%d.txt",carpeta_output, ult + i);
        evolucion_hasta_decir_basta(filename_input, N_pasos, dt, K, betta,filename_output);
    }

    // Guardar historial
    char historial_path[512];
    sprintf(historial_path, "%s\\ER_Historial.txt", carpeta_output);  

    FILE* historial = fopen(historial_path, "a");
    if (historial == NULL) {
        perror("No se pudo abrir el archivo de historial");
        return;
    }

    time_t t = time(NULL);
    struct tm tm = *localtime(&t);

    fprintf(historial, "--------------------------------------------------------------------------------\n");
    fprintf(historial, "Fecha: %04d-%02d-%02d %02d:%02d:%02d\n", 
            tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday,
            tm.tm_hour, tm.tm_min, tm.tm_sec);
    fprintf(historial, "Redes (input): ER_%d.txt a ER_%d.txt\n", rede_ini, rede_ini + N_redes - 1);
    fprintf(historial, "Resultados (output): %s\n", filename_output);
    fprintf(historial, "Parámetros -> K: %.2lf | β: %.2lf | dt: %.4lf | N_pasos: %.0lf\n", 
            K, betta, dt, N_pasos);
    fprintf(historial, "--------------------------------------------------------------------------------\n\n");

    fclose(historial);
}



void evolucion_hasta_decir_basta(char*filename_input, int N_pasos, double dt, double K, double betta,char*filename_output){
    int* vecinos;
    int* grados;
    int total_nodos = 0;
    double delta;
    double op_media, desvest;
    leer_red(filename_input, &vecinos, &grados, &total_nodos);
    double* x = malloc(total_nodos * sizeof(double));
    delta= condiciones_iniciales(K,total_nodos,x);
    bool flag=false;
    FILE *archivo;
    archivo = fopen(filename_output, "w"); //Para crear el archivo
    fprintf(archivo, "%lf\t%lf\t%d\t%lf\t%s\n",K, betta, N_pasos, dt,filename_input);
    if (archivo == NULL) {
        fprintf(stderr, "Error al abrir el archivo %s\n", filename_output);
        return;
    }
    fclose(archivo);

    while(flag==false){
    for(int j=0;j<N_pasos;j++){
        rk4_step(x,total_nodos,dt,K, betta,delta,vecinos, grados);
        j++;
    }
    flag= esta_termalizada(K,betta,delta,vecinos,grados,x,total_nodos);
    }
    op_media=valor_medio(x,total_nodos);
    desvest=desviacion_estandar(x,total_nodos,op_media);

    archivo = fopen(filename_output, "a");  // "a" para añadir sin sobrescribir
    fprintf(archivo, "%lf\t%lf\n",op_media,desvest);
    for(int i=0;i<total_nodos;i++){
        fprintf(archivo, "%lf\n", x[i]);
    }
    fclose(archivo);
}

void calcular_fraccion_polarizados(double K, double betta, int N_res, const char* filename_output) {
    char carpeta_path[128];
    char filepath[256];
    FILE* archivo;
    int polarizados = 0;
    int no_polarizados = 0;

    // Construir ruta de carpeta: carpeta K\betta
    sprintf(carpeta_path, "%.2lf\\%.2lf", K, betta);

    for (int i = 0; i < N_res; i++) {
        // Construir ruta completa del archivo ER_i.txt
        sprintf(filepath, "EJECUTABLES\\%s\\ER_%d.txt", carpeta_path, i);
        archivo = fopen(filepath, "r");
        if (archivo == NULL) {
            fprintf(stderr, "⚠️ No se pudo abrir el archivo: %s\n", filepath);
            continue;
        }

        char linea[512];
        // Leer y descartar la primera línea
        if (fgets(linea, sizeof(linea), archivo) == NULL) {
            fclose(archivo);
            continue;
        }

        // Leer la segunda línea
        // Leer la segunda línea
        if (fgets(linea, sizeof(linea), archivo) != NULL) {
            double media, desvest;
         if (sscanf(linea, "%lf %lf", &media, &desvest) == 2) {
            
              if(fabs(media)<0.5 && desvest<0.5){
                 no_polarizados++;
              }
               else{
                    if (desvest>media){
                       polarizados++;
                  }
                 else{
                    no_polarizados++;
                }
            }
        

        }
    }


        fclose(archivo);
    }

    // Calcular fracciones
    double frac_polarizados = (double)polarizados / N_res;
    double frac_no_polarizados = (double)no_polarizados / N_res;

    // Escribir en archivo de salida
    FILE* out = fopen(filename_output, "a");
    if (out == NULL) {
        fprintf(stderr, "❌ No se pudo abrir el archivo de salida: %s\n", filename_output);
        return;
    }

    fprintf(out, "%.2lf %.2lf %.6lf %.6lf\n", K, betta, frac_polarizados, frac_no_polarizados);
    fclose(out);
}
