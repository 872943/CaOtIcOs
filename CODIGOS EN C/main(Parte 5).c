#include "funciones.h"
/*  Notad que, si el grado es lo suficientemente pequeño (y por tanto el horizonte
    de información), puede aparecer polarización incluso en ausencia de homofilia
    (beta = 0). Estudiad sistemáticamente esta aparición con redes ER de grado bajo,
    y con redes WS de grado k= 10 con diversas probabilidades de rewiring
    (reproduciendo cualitativamente la figura 4 del paper, pero con fracciones de
    polarización).*/    


int main() {
    
    
        srand(time(NULL)); // Inicializa aleatoriedad
    
        double K=10;         // a partir de ahora fijamos K=10
        double betta=0;

        int N_simulaciones = 100;
        int N_pasos=280;   
        double dt=0.1;

       