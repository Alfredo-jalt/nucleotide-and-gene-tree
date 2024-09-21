#José Alfredo López Torres 

#Reto con viruses de diferentes países

library("seqinr")
#Sacamos la carpeta de base de datos
getwd()
setwd("Workspace url")
#Leemos todas las variables de genomas
wuhanhu <- read.fasta("URL")
argentina <- read.fasta("URL")
australia <- read.fasta("URL")
brazil <- read.fasta("URL")
canada <- read.fasta("URL")
colombia <- read.fasta("URL")
france <- read.fasta("URL")
germany <- read.fasta("URL")
guatemala <- read.fasta("URL")
india <- read.fasta("URL")
italy <- read.fasta("URL")
japan <- read.fasta("URL")
korea <- read.fasta("URL")
mexico <- read.fasta("URL")
netherlands <- read.fasta("URL")
south_africa <- read.fasta("URL")
spain <- read.fasta("URL")
taiwan <- read.fasta("URL")
turkey <- read.fasta("URL")
usa <- read.fasta("URL")


# List of variable names and associated data
viruses <- list(
  wuhanhu = wuhanhu,
  argentina = argentina,
  australia = australia,
  brazil = brazil,
  canada = canada,
  colombia = colombia,
  france = france,
  germany = germany,
  guatemala = guatemala,
  india = india,
  italy = italy,
  japan = japan,
  korea = korea,
  mexico = mexico,
  netherlands = netherlands,
  south_africa = south_africa,
  spain = spain,
  taiwan = taiwan,
  turkey = turkey,
  usa = usa
)

#Funcion para sacar la longitud de los viruses
print_virus_lengths <- function(viruses) {
  for (virus_name in names(viruses)) {
    cat(virus_name, "\n")
    sequences <- viruses[[virus_name]]  # Acceder a las secuencias dentro del objeto del virus
    total_length <- length(sequences[[1]])  # Sacar la longitud de las secuencias
    cat("Total length of the genome: ", total_length, "\n\n")
  }
}

# Llamada a la función 
print_virus_lengths(viruses)

######################################################


# Instala y carga la librería necesaria
#install.packages("ggplot2")
library(ggplot2)

# Función para calcular la frecuencia de caracteres en los genomas de virus
calculate_nucleotide_frequency <- function(viruses) {
  frequencies <- list()
  for (virus_name in names(viruses)) {
    sequences <- viruses[[virus_name]]
    nucleotides <- c("a", "g", "c", "t")
    counts <- sapply(nucleotides, function(nucleotide) sum(sapply(sequences, function(sequence) sum(sequence == nucleotide))))
    frequencies[[virus_name]] <- counts
  }
  return(frequencies)
}

# Llama a la función para calcular las frecuencias
nucleotide_frequencies <- calculate_nucleotide_frequency(viruses)

# Organiza los datos para crear el dataframe
country_names <- names(nucleotide_frequencies)
nucleotides <- c("a", "g", "c", "t")
data <- data.frame()
for (country in country_names) {
  country_data <- data.frame(country = country, nucleotide = nucleotides, frequency = nucleotide_frequencies[[country]])
  data <- rbind(data, country_data)
}

# Crea la gráfica
ggplot(data, aes(x = country, y = frequency, fill = nucleotide)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Frecuencia de nucleótidos por país",
       x = "País",
       y = "Frecuencia") +
  scale_fill_manual(values = c("a" = "red", "g" = "blue", "c" = "green", "t" = "purple")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##############################################################################

library("Rcpp")
library("ape")
library("phangorn")
library("phytools")
library("geiger")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")
library("Biostrings")
library("seqinr")
library("adegenet")
library("ape")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")
library("ggtree")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DECIPHER")
library("DECIPHER")
library("viridis")
library("ggplot2")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggmsa")
library("ggmsa")

# Definir nombres personalizados para las secuencias
nombres <- c("Argentina", "Australia", "Brasil", "Canada", "Colombia", "France",
             "Germany", "Guatemala", "India", "Italy", "Japan", "Korea",
             "Mexico", "Netherlands", "South Africa", "Spain", "Taiwan", "Turkey",
             "USA", "WuhanHu China")

virus <- c("OP302820.1", "OP848496.1", "OP855522.1", "OR999078.1", "OQ551269.1", "OQ890287.1",
           "OQ503465.1", "OP314044.1", "OQ852631.1", "OR761953.1", "OQ804217.1", "OR447552.1",
           "OR157741.1", "OR427989.1", "OR939737.1", "ON115272.1", "OQ931803.1", "OR529199.1",
           "PP064475.1", "NC_045512.2")

secuencias <- read.GenBank(virus)

# Asignar nombres personalizados
names(secuencias) <- nombres

# Crear fasta de las secuencias juntas
write.dna(secuencias,  file = "virus_seqs.fasta", format = "fasta", append = FALSE, nbcol = 6, colsep = " ", colw = 10)

virus_not_align <- readDNAStringSet("virus_seqs.fasta", format = "fasta")
virus_not_align

# Organizar secuencias
virus_not_align <- OrientNucleotides(virus_not_align)

# Alineamiento
virus_seq_align <- AlignSeqs(virus_not_align)

# Guardar resultado
writeXStringSet(virus_seq_align, file = "virus_seq_align.fasta")

# Leer nuevo Fasta
virus_aligned <- read.alignment("virus_seq_align.fasta", format = "fasta")

# Matriz de distancia
matriz_distancia <- dist.alignment(virus_aligned, matrix = "similarity")

temp <- as.data.frame(as.matrix(matriz_distancia))
table.paint(temp, cleg = 0, clabel.row = 0.5, clabel.col = 0.5) + scale_color_viridis()

virus_tree <- nj(matriz_distancia)
class(virus_tree) 

virus_tree <- ladderize(virus_tree)

# Generar Arbol
plot(virus_tree, cex = 0.6)
title("Arbol de los diferentes virus")

# Plot utilizando ggtree que es parte de ggplot:
ggtree(virus_tree) + geom_tiplab()
ggtree(virus_tree, branch.length = 'none', layout = 'circular') + geom_tiplab()

