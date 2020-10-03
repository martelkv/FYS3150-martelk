# Prosjekt 2 - FYS3150
Navn: Marte Lunde Kvam

## Hvordan kjøre programmene
1. Hvis du har "make" installert, kjør ```make``` i terminalen når du er i mappen. Dette vil kompilere filen.  

2. For å kjøre kildekoden til oppgave 2b, kjør ```make run_jacobi```. For å få figur 3 og 4 fra rapporten, endre verdien i run_jacobi makefile til 10 (for figur 3) til 100 for figur 4, og kjør ```make run_jacobi```. Vi oppretter også da fil *plott_egenvektor.txt*, hvor egenvektorene lagres. *plott_egenvektor.py* vil også kjøres når man kjører kommandoen over.  

3. For å kjøre Armadillo-versjonen av oppgave 2b, kjør ```make run_armendillo```.  I *makefile* kan man endre antall integrasjonspunkter ved å endre argumentet som sendes inn med *./armendillo_eig.exe*, f. eks. 1000. Dette argumentet må være med for at programmet kan kjøres. For resultatene i oppgaven, kjørte jeg med 10 integrasjonspunkter.  

4. For resultatene i oppgave 2d og 2e, kjør ```make run_elektroner```. I *makefile* kan man endre antall integrasjonspunkter ved å endre argumentet som sendes inn med *./quantum_dots.exe*, f. eks. 1000. Dette argumentet må være med for at programmet kan kjøres. For resultatene i oppgaven, kjørte jeg med 100 integrasjonspunkter.    

## Avhengigheter
1. For å kompilere programmet, burde man ha *make* installert. Hvis man ikke har dette installert, kan man skrive inn kommandoene i *makefile* etter hverandre.  
