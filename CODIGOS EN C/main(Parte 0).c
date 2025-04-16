#include "funciones.h"


int main() {
    srand(time(NULL)); // Inicializa aleatoriedad
    int N_sim=100;
    int N_pasos=280;
    double dt=0.1;
    double K=10;
    double betta=0.1;
    int number_name=99;

    muchas_simulaciones_ER(N_sim, N_pasos, dt,K, betta, number_name);
}