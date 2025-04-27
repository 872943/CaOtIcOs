#ifndef FUNCIONES_H
#define FUNCIONES_H

// Librerías necesarias
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <dirent.h>
#include <ctype.h>
#include <stdbool.h>
#include <unistd.h> // Para sleep en Linux, si es necesario
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#ifdef _WIN32
#include <direct.h>
#define mkdir_p(dir) _mkdir(dir)
#else
#include <unistd.h>
#define mkdir_p(dir) mkdir(dir, 0777)
#endif

// Declaraciones de funciones
void leer_red(const char* filename, int** vecinos, int** grados, int* total_nodos);

double sumando_i(double K, double betta, double delta, int* vecinos, int grado, double* x, int i);

void derivada(double K, double betta, double delta, int* vecinos, int *grados, double* x, double* derivada_x, int total_nodos);

void rk4_step(double* x, int total_nodos, double dt, double K, double betta, double delta, int* vecinos, int* grados);

double condiciones_iniciales(double K, int total_nodos, double *x0);

void polarizacion(double *x, int total_nodos, double *op_media, double *desvest);

void copia_conf(double* x_original, double *x_copiado, int total_nodos);

void escribe_evolucion_polarizacion(double tiempo, double x_med, double x_desvest, const char* filename_output);

void simula_polarizacion(const char* filename_input, int N_pasos, double dt, double K, double betta, const char* filename_output);

int obtener_siguiente_indice(const char* carpeta);

void muchas_simulaciones_ER(int N_sim, int N_pasos, double dt, double K, double betta, int number_name);

void muchas_simulaciones_WS(int N_sim, int N_pasos, double dt, double K, double betta, int number_name);

bool es_polarizada(const char* filename);

double valor_medio(double arr[], int tam);

double desviacion_estandar(double arr[], int tam, double media);

void escribe_evolucion_individual(char* filename_output, int N_nodos, double *x, double tiempo);

void escribe_velocidad_modulo(char* filename_output, double tiempo, double velocidad);

void evolucion_persona_a_persona(char* filename_input, char* filename_output, int N_pasos, double dt, double K, double betta, char* filename_output_velocidad);

double cacula_velocidad_modulo(double K, double betta, double delta, int* vecinos, int *grados, double* x, int total_nodos);

bool esta_termalizada(double K, double betta, double delta, int* vecinos, int* grados, double* x, int total_nodos);

void frac_polarizado(int N_redes, int rede_ini, double K, double betta, double N_pasos, double dt);

void evolucion_hasta_decir_basta(char*filename_input, int N_pasos, double dt, double K, double betta,char*filename_output);

void calcular_fraccion_polarizados(double K, double betta, int N_res, const char* filename_output);

void evolucion_hasta_decir_basta_apartado5(char*filename_input, int N_pasos, double dt, double k,char*filename_output);

void frac_polarizado_apartado5 (int N_redes, int rede_ini, double k, double N_pasos, double dt);

void calcular_fraccion_polarizados_apartado5(double k, int N_res, const char* filename_output);

double leer_k_desde_parametros(int red_id) ;
#endif
