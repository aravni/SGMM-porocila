# SGMM-porocila
Ta repozitorij vsebuje poročili za vaje pri predmetu Sestavljena gonila v mobilni tehniki, Fakulteta za strojništvo Univerze v ljubljani, letnik 2024/25.

## Smotrna raba programske kode
Koda, spisana v programskem jeziku Python, je namenjena v pomoč pri pisanju poročila. Nikakor ne želim vzpodbujati prepisovanja ali uporabe kode za generiranje slik. Uporabnik je vljudno povabljen tako k nadgradnji, kot tudi iskanju napak v kodi in dodajanju podatkov za nova vozila.

## Sestava repozitorija

1) Direktorij `podatki/` z datotekami:
- Primeri datotek s podatki o vozilu. Za analizo novega vozila je potrebno spremeniti podatke v datoteki, vendar ohraniti strukturo le-te.

2) Programske datoteke za 1. poročilo:
- `porocilo1.py` za obdelavo podatkov vozila iz datoteke, tvorijo se osnovni grafi brez upoštevanja zdrsa.
- `ReadData.py` vsebuje funkcijo za branje podatkov vozila iz datoteke.
- `porocilo1_zdrs.py` obravnava poleg bistvenih zadev in kode `porocilo1.py` tudi silo trenja in zdrs pogonskih koles. Na tem mestu so že možne nadgradnje.
- `porocilo1.opti.py` poizkuša prilagoditi proste parametre simulacije ali prestavna razmerja, da je dosežen željen čas do 100 km/h in željena maksimalna hitrost. Tukaj je še največ prostora za nadgranjo, saj je ta koda še v delu in ni dodelana do željene mere za oddajo.

3) Pomožne datoteke:
- README.md, LICENSE

## Opombe
- Koda je preverjena in deluje v `Python 3.11.5`
- V kolikor uporabnik nima dodane custom kode za poravnavo in stil grafov, je potrebno le-te ali naložiti kot module, ali opustiti v programski kodi, vendar v slednjem primeru pravilen izris grafov ni zagotovljen.
