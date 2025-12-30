# Approssimazione di un Sistema Differenziale Non Lineare

Questo progetto analizza le soluzioni di un sistema di equazioni differenziali del primo ordine non lineare tramite diversi metodi numerici.
Il sistema studiato è tridimensionale e vincolato a una superficie sferica. L'obiettivo è analizzare le traiettorie e stimare l'errore di approssimazione partendo da dati iniziali sulla sfera unitaria.

## Metodi Implementati
Il software include l'implementazione in C++ dei seguenti metodi di ordine $p=4$:
* **Gauss 2**: Metodo implicito isometrico che conserva la norma.
* **Adam-Bashforth 4**: Metodo multistep esplicito.
* **Runge-Kutta 4**: Metodo esplicito ad un passo.
e l'implementazione del seguente metodo di ordine $p=6$:
* **Gauss 3**: Metodo implicito isometrico.

