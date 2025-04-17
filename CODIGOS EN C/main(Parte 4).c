#include "funciones.h"


int main() {
//   srand(time(NULL)); // Inicializa aleatoriedad
//   int N_sim=100;
//   int N_pasos=280;
//   double dt=0.1;
//   double betta=0.1;
//   int number_name=99;

//   muchas_simulaciones_ER(N_sim, N_pasos, dt,K, betta, number_name);


    srand(time(NULL)); // Inicializa aleatoriedad

    double K=10;         // a partir de ahora fijamos K=10
    int N_sim=200;
    int N_pasos=280;   
    double dt=0.1;
    double betta_max=1.0
    double betta=0
    double N, delta_betta
    int N_sim=200
    N=100
    delta_betta=(betta_max -betta)/N
    //empezamos con B=0 y realizamos los caluculos, aumentamos la betta un
    // delta_betta y reprtimos el proceso, y asi sucesivamente hasta B=1.

    //queremos ver cuál es la fracción de estados en los que obtenemos polarización.
    //Para ello, vamos a simular 200 veces para cada valor de betta para cada <k>.
    //las redes ER con <k>=6 son las que van de ER_200 a ER_399
    //las redes ER con <k>=14 son las que van de ER_400 a ER_599

    int number_name=200;
    
    //bucle para <k>=6
    for(int i=0,i<N,i++) {
        muchas_simulaciones_ER(N_sim, N_pasos, dt,K, betta, number_name);
        betta=betta+delta_betta
    }


    int number_name=400;
    betta=0

    //bucle para <k>=14
    for(int i=0,i<N,i++) {
        muchas_simulaciones_ER(N_sim, N_pasos, dt,K, betta, number_name);
        betta=betta+delta_betta
    }




}