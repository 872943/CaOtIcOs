#include "funciones.h"


int main() {

    srand(time(NULL)); // Inicializa aleatoriedad
    int N_sim=100;
    int N_pasos=1000;
    double dt=0.1;
    double K=10;
    double betta=1.5;
    int number_name=99;
    int total_nodos;

    int* grados = NULL;
    int* vecinos = NULL;



       
    evolucion_persona_a_persona(
        "ARCHIVOS_REDES/ER/ER_0.txt",
        "Resultados (PARTE 0)/Evolucion individual/ER/evoluciones/ER_0.txt",
        N_pasos,
        dt,
        K,
        betta,
        "Resultados (PARTE 0)/Evolucion individual/ER/velocidades/ER_0.txt"
    );

    return 0;
}