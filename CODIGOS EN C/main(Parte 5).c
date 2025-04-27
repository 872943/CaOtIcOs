
/*  Notad que, si el grado es lo suficientemente pequeño (y por tanto el horizonte
    de información), puede aparecer polarización incluso en ausencia de homofilia
    (beta = 0). Estudiad sistemáticamente esta aparición con redes ER de grado bajo,
    y con redes WS de grado k= 10 con diversas probabilidades de rewiring
    (reproduciendo cualitativamente la figura 4 del paper, pero con fracciones de
    polarización).*/    

#include "funciones.h"


int main() {
    srand(time(NULL)); // Inicializa aleatoriedad
    double dt = 0.1;
    int N_pasos = 50;
    int N_redes = 1;
    int rede_ini = 745;    // Comenzar desde ER_745
    int rede_final = 805;  // Terminar en ER_805
    int N_sim = 100;
    double delta_k = 0.1;
    double k_actual = 2.0;
    double k_final = 8.0;

    int id_red = rede_ini;

    while (id_red <= rede_final && k_actual <= k_final) {
        // Leer el valor real de k del archivo de parámetros
        double k_real = leer_k_desde_parametros(id_red);

        printf("Procesando red ER_%d con k_real = %.4lf\n", id_red, k_real);

        // Llamadas a las funciones usando k_real
        frac_polarizado_apartado5(N_redes, id_red, k_real, N_pasos, dt);
        calcular_fraccion_polarizados_apartado5(k_real, N_redes, "Resultados.txt");

        // Avanzar
        id_red++;
        k_actual += delta_k;
    }

    return 0;
}