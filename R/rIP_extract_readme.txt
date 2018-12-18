GUIDA a rIP extract.
tutti gli esempi usano la mediana come FUN, PMBest-sync come stream, e red/blue/green come livelli della categoria "color"

epochStreamApply 
per ciascuna istanza di ciascun colore calcola la mediana.
salva i risultati nel DyadCategory session$color

catExtract(experiment, category="color", by, FUN = median, ...)
prende tutte le tabelle categorie di tutte le sedute, le unisce e (se specificato) le aggrega con la mediana (ad es per seduta, o per colore)
Funziona in combo con  epochStreamApply calcolando ad esempio la mediana di tutte le mediane di "red" delle diverse occorrenze


epochStreamApplyLong
unisce ciascuna istanza di ciascun colore in un unico stream.
salva i risultati nel DyadSignal: session$signal$color_PMBestSync_Median -> list("red"=... , "blue"=..., "green"=...)


streamExtract
prende tutti gli stream PMBest-sync e applica la mediana


################### TODO
in generale quando voglio un summary delle mie analisi SYNC voglio:
- estrarre tutto o tutte le finestre corrispondenti ad una certa categoria, in un certo stream.
- calcolarne una funzione riassuntiva( media/mediana) per finestra o incollando tutte le finestre fra loro (per sessione o per esperimento)
- estrarre i valori in una tabella, aggregando secondo diversi criteri

summary.Dyad... <- function(x,
						epoch = c({missing}, "IM", "PACS"),
						stream=c("PMBest","sync"),
						aggregateBy=c("session", "dyad",...), #passed to aggregate
						epochBy = c("each","session","dyad","experiment"), #epoche per istanza o incollate?
						summaryFUN = "mean", #che funzione usare per aggregare
						epochFUN = "median" #con che funzione estrarre 1 valore da ogni singola finestra (o collezione di finestre)
						)
						
