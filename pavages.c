/*
 * pavages — résoudre des puzzles de pavage avec l'algorithme X de Knuth
 * et les « liens dansants » (Dancing Links, DLX).
 *
 * Un puzzle de pavage est un problème de COUVERTURE EXACTE : une matrice de
 * 0 et de 1, et il faut choisir des lignes qui ont exactement un 1 dans
 * chaque colonne. Ici :
 *   - une colonne par pièce (chaque pièce sert une fois) et une par case ;
 *   - une ligne par façon de poser une pièce (orientation + position).
 *
 * Puzzles fournis : les 12 pentaminos dans un rectangle (6x10, 5x12, 4x15,
 * 3x20), dans un carré 8x8 troué, dans une forme lue dans un fichier, et le
 * cube Soma en 3D (7 pièces, cube 3x3x3).
 *
 * Les solutions qui ne diffèrent que par une symétrie du plateau (rotation,
 * miroir) sont regroupées : chaque solution est ramenée à une forme
 * canonique et rangée dans une table de hachage.
 */
#define _POSIX_C_SOURCE 200809L
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

static void meurt(const char *m) {
    fprintf(stderr, "pavages : %s\n", m);
    exit(1);
}
static void *alloue(size_t n) {
    void *p = calloc(1, n ? n : 1);
    if (!p) meurt("plus de mémoire");
    return p;
}
static void *agrandit(void *p, size_t n) {
    p = realloc(p, n);
    if (!p) meurt("plus de mémoire");
    return p;
}

/* ============================================================== DLX */

typedef struct { int g, d, h, b, col, ligne; } Noeud;

static Noeud *N;
static int nb_noeuds, cap_noeuds, nb_cols, *taille;
static long long mises_a_jour; /* « updates » au sens de Knuth : liens défaits */

static int nouveau_noeud(void) {
    if (nb_noeuds == cap_noeuds) {
        cap_noeuds = cap_noeuds ? cap_noeuds * 2 : 4096;
        N = agrandit(N, sizeof(Noeud) * cap_noeuds);
    }
    return nb_noeuds++;
}

/* le nœud 0 est la racine, les nœuds 1..n sont les en-têtes de colonnes */
static void dlx_init(int ncols) {
    nb_noeuds = 0;
    nb_cols = ncols;
    taille = alloue(sizeof(int) * (ncols + 1));
    for (int i = 0; i <= ncols; i++) {
        int k = nouveau_noeud();
        N[k] = (Noeud){(k + ncols) % (ncols + 1), (k + 1) % (ncols + 1), k, k, k, -1};
    }
}

static void dlx_ligne(const int *cols, int k, int id) {
    int premier = -1;
    for (int i = 0; i < k; i++) {
        int c = cols[i] + 1, x = nouveau_noeud();
        N[x].col = c;
        N[x].ligne = id;
        N[x].h = N[c].h; N[x].b = c;
        N[N[c].h].b = x; N[c].h = x;
        taille[c]++;
        if (premier < 0) { premier = x; N[x].g = N[x].d = x; }
        else { N[x].d = premier; N[x].g = N[premier].g; N[N[premier].g].d = x; N[premier].g = x; }
    }
}

static void couvre(int c) {
    N[N[c].d].g = N[c].g;
    N[N[c].g].d = N[c].d;
    for (int i = N[c].b; i != c; i = N[i].b)
        for (int j = N[i].d; j != i; j = N[j].d) {
            N[N[j].b].h = N[j].h;
            N[N[j].h].b = N[j].b;
            taille[N[j].col]--;
            mises_a_jour++;
        }
}

static void decouvre(int c) {
    for (int i = N[c].h; i != c; i = N[i].h)
        for (int j = N[i].g; j != i; j = N[j].g) {
            taille[N[j].col]++;
            N[N[j].b].h = j;
            N[N[j].h].b = j;
        }
    N[N[c].d].g = c;
    N[N[c].g].d = c;
}

typedef void (*Rappel)(const int *lignes, int k);
static int choix[64];

static void cherche(int k, Rappel rappel) {
    if (N[0].d == 0) { rappel(choix, k); return; }
    /* heuristique de Knuth : la colonne la moins remplie (moins de branches) */
    int c = N[0].d;
    for (int j = N[c].d; j != 0; j = N[j].d)
        if (taille[j] < taille[c]) c = j;
    if (taille[c] == 0) return;
    couvre(c);
    for (int r = N[c].b; r != c; r = N[r].b) {
        choix[k] = N[r].ligne;
        for (int j = N[r].d; j != r; j = N[j].d) couvre(N[j].col);
        cherche(k + 1, rappel);
        for (int j = N[r].g; j != r; j = N[j].g) decouvre(N[j].col);
    }
    decouvre(c);
}

/* ======================================================= géométrie 3D */

typedef struct { int x, y, z; } P3;
typedef struct { int m[3][3]; int det; } Mat;

static Mat mats[48];
static int nb_mats;

/* les 48 matrices de permutation signées : toutes les isométries du cube */
static void genere_mats(void) {
    static const int perms[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    static const int signe_perm[6] = {1, -1, -1, 1, 1, -1};
    for (int p = 0; p < 6; p++)
        for (int s = 0; s < 8; s++) {
            Mat M = {{{0}}, 0};
            int det = signe_perm[p];
            for (int i = 0; i < 3; i++) {
                int sg = (s >> i & 1) ? -1 : 1;
                M.m[i][perms[p][i]] = sg;
                det *= sg;
            }
            M.det = det;
            mats[nb_mats++] = M;
        }
}

static P3 applique(const Mat *M, P3 p) {
    int v[3] = {p.x, p.y, p.z}, r[3];
    for (int i = 0; i < 3; i++) r[i] = M->m[i][0] * v[0] + M->m[i][1] * v[1] + M->m[i][2] * v[2];
    return (P3){r[0], r[1], r[2]};
}

static int cmp_p3(const void *a, const void *b) {
    const P3 *p = a, *q = b;
    if (p->z != q->z) return p->z - q->z;
    if (p->y != q->y) return p->y - q->y;
    return p->x - q->x;
}

/* translate pour que le minimum soit en (0,0,0), puis trie : forme normale */
static void normalise(P3 *c, int n) {
    int mx = c[0].x, my = c[0].y, mz = c[0].z;
    for (int i = 1; i < n; i++) {
        if (c[i].x < mx) mx = c[i].x;
        if (c[i].y < my) my = c[i].y;
        if (c[i].z < mz) mz = c[i].z;
    }
    for (int i = 0; i < n; i++) { c[i].x -= mx; c[i].y -= my; c[i].z -= mz; }
    qsort(c, n, sizeof(P3), cmp_p3);
}

/* ================================================================ puzzle */

#define MAX_CUBES 6
#define MAX_ORIENT 48
typedef struct {
    char nom;
    int n;
    P3 cubes[MAX_CUBES];
    int nb_orient;
    P3 orient[MAX_ORIENT][MAX_CUBES];
} Piece;

typedef struct {
    int L, H, P;        /* dimensions de la boîte englobante */
    int *num;           /* num[(z*H + y)*L + x] = numéro de case, ou -1 */
    P3 *cases;
    int nb_cases;
    Piece pieces[16];
    int nb_pieces;
    int plan;           /* 1 : puzzle plat (pentaminos), 0 : 3D (Soma) */
    char miroir_a, miroir_b; /* pièces échangées par un miroir (Soma : A et B) */
    /* symétries du plateau : permutations des cases (+ miroir ?) */
    int nb_sym;
    int *sym[48];
    int sym_miroir[48];
} Puzzle;

static Puzzle Z;

static int case_en(int x, int y, int z) {
    if (x < 0 || y < 0 || z < 0 || x >= Z.L || y >= Z.H || z >= Z.P) return -1;
    return Z.num[(z * Z.H + y) * Z.L + x];
}

static void plateau_depuis_texte(const char *texte, int profondeur) {
    /* chaque ligne est une rangée ; '.' ou ' ' = trou ; couches séparées par une ligne vide */
    int L = 0, H = 0, x = 0;
    for (const char *p = texte; *p; p++) {
        if (*p == '\n') { if (x > L) L = x; x = 0; H++; }
        else x++;
    }
    if (x) { if (x > L) L = x; H++; }
    H /= profondeur;
    Z.L = L; Z.H = H; Z.P = profondeur;
    Z.num = alloue(sizeof(int) * L * H * profondeur);
    Z.cases = alloue(sizeof(P3) * L * H * profondeur);
    for (int i = 0; i < L * H * profondeur; i++) Z.num[i] = -1;
    int y = 0;
    x = 0;
    for (const char *p = texte; *p; p++) {
        if (*p == '\n') { x = 0; y++; continue; }
        if (*p != '.' && *p != ' ' && *p != '\r' && y < H * profondeur) {
            int z = y / H, yy = y % H;
            Z.num[(z * H + yy) * L + x] = Z.nb_cases;
            Z.cases[Z.nb_cases++] = (P3){x, yy, z};
        }
        x++;
    }
}

static void ajoute_piece(char nom, const char *forme) {
    /* forme : rangées séparées par '/', '#' = cube ; couches séparées par '|' */
    Piece *p = &Z.pieces[Z.nb_pieces++];
    p->nom = nom;
    int x = 0, y = 0, z = 0;
    for (const char *c = forme; *c; c++) {
        if (*c == '/') { y++; x = 0; }
        else if (*c == '|') { z++; y = 0; x = 0; }
        else { if (*c == '#') p->cubes[p->n++] = (P3){x, y, z}; x++; }
    }
}

/* toutes les orientations distinctes : 8 dans le plan, ou les 24 rotations en 3D */
static void orientations(Piece *p) {
    for (int m = 0; m < nb_mats; m++) {
        const Mat *M = &mats[m];
        if (Z.plan ? M->m[2][2] != 1 : M->det != 1) continue;
        P3 t[MAX_CUBES];
        for (int i = 0; i < p->n; i++) t[i] = applique(M, p->cubes[i]);
        normalise(t, p->n);
        int deja = 0;
        for (int o = 0; o < p->nb_orient && !deja; o++) deja = !memcmp(p->orient[o], t, sizeof(P3) * p->n);
        if (!deja) memcpy(p->orient[p->nb_orient++], t, sizeof(P3) * p->n);
    }
}

/* symétries du plateau : celles des 48 isométries qui l'envoient sur lui-même */
static void symetries(void) {
    for (int m = 0; m < nb_mats; m++) {
        const Mat *M = &mats[m];
        if (Z.plan && M->m[2][2] != 1) continue;
        int *perm = alloue(sizeof(int) * Z.nb_cases), ok = 1;
        for (int c = 0; c < Z.nb_cases && ok; c++) {
            /* coordonnées doublées et centrées, pour que le centre soit entier */
            P3 q = {2 * Z.cases[c].x - (Z.L - 1), 2 * Z.cases[c].y - (Z.H - 1), 2 * Z.cases[c].z - (Z.P - 1)};
            P3 r = applique(M, q);
            int x2 = r.x + (Z.L - 1), y2 = r.y + (Z.H - 1), z2 = r.z + (Z.P - 1);
            if ((x2 | y2 | z2) & 1) { ok = 0; break; }
            perm[c] = case_en(x2 / 2, y2 / 2, z2 / 2);
            if (perm[c] < 0) ok = 0;
        }
        if (ok) {
            Z.sym[Z.nb_sym] = perm;
            Z.sym_miroir[Z.nb_sym++] = M->det < 0 && !Z.plan;
        } else free(perm);
    }
}

/* ============================================= table des placements */

typedef struct { int piece; int cases[MAX_CUBES]; } Placement;
static Placement *placements;
static int nb_placements;

/* garde un placement seulement s'il est le plus petit de son orbite sous les
 * symétries : supprime les solutions symétriques sans en perdre aucune classe */
static int minimal_dans_orbite(const int *cases, int n) {
    int tri[MAX_CUBES], img[MAX_CUBES];
    memcpy(tri, cases, sizeof(int) * n);
    for (int i = 1; i < n; i++)
        for (int j = i; j > 0 && tri[j] < tri[j - 1]; j--) { int t = tri[j]; tri[j] = tri[j - 1]; tri[j - 1] = t; }
    for (int s = 0; s < Z.nb_sym; s++) {
        for (int i = 0; i < n; i++) img[i] = Z.sym[s][tri[i]];
        for (int i = 1; i < n; i++)
            for (int j = i; j > 0 && img[j] < img[j - 1]; j--) { int t = img[j]; img[j] = img[j - 1]; img[j - 1] = t; }
        for (int i = 0; i < n; i++) {
            if (img[i] < tri[i]) return 0;
            if (img[i] > tri[i]) break;
        }
    }
    return 1;
}

static void construit_matrice(int briser_symetrie) {
    int cap = 1024;
    placements = alloue(sizeof(Placement) * cap);
    dlx_init(Z.nb_pieces + Z.nb_cases);
    /* la pièce qui sert à briser la symétrie : celle qui a le moins d'orientations */
    int casseuse = -1;
    if (briser_symetrie) {
        for (int i = 0; i < Z.nb_pieces; i++)
            if (casseuse < 0 || Z.pieces[i].nb_orient < Z.pieces[casseuse].nb_orient) casseuse = i;
        /* il faut une pièce entièrement symétrique et un jeu sans paire chirale */
        if (Z.pieces[casseuse].nb_orient != 1 || Z.miroir_a) casseuse = -1;
    }
    for (int i = 0; i < Z.nb_pieces; i++) {
        Piece *p = &Z.pieces[i];
        for (int o = 0; o < p->nb_orient; o++)
            for (int dz = 0; dz < Z.P; dz++)
                for (int dy = 0; dy < Z.H; dy++)
                    for (int dx = 0; dx < Z.L; dx++) {
                        int cols[MAX_CUBES + 1], ok = 1;
                        Placement pl = {i, {0}};
                        cols[0] = i;
                        for (int k = 0; k < p->n && ok; k++) {
                            int c = case_en(p->orient[o][k].x + dx, p->orient[o][k].y + dy, p->orient[o][k].z + dz);
                            if (c < 0) ok = 0;
                            else { pl.cases[k] = c; cols[k + 1] = Z.nb_pieces + c; }
                        }
                        if (!ok) continue;
                        if (i == casseuse && !minimal_dans_orbite(pl.cases, p->n)) continue;
                        if (nb_placements == cap) placements = agrandit(placements, sizeof(Placement) * (cap *= 2));
                        placements[nb_placements] = pl;
                        dlx_ligne(cols, p->n + 1, nb_placements++);
                    }
    }
}

/* ======================================= solutions et forme canonique */

static char *ensemble;      /* table de hachage de chaînes de taille nb_cases */
static size_t cap_ens, nb_ens;
static long long nb_brutes, nb_total;
static int a_afficher, couleurs;
static char *premiere;
static char *gardees[64];
static int nb_gardees;

static uint64_t hache(const char *s, int n) {
    uint64_t h = 1469598103934665603ull;
    for (int i = 0; i < n; i++) h = (h ^ (uint8_t)s[i]) * 1099511628211ull;
    return h;
}

static int insere(const char *s) {
    int n = Z.nb_cases;
    if (2 * (nb_ens + 1) > cap_ens) {
        size_t ancien = cap_ens;
        char *vieux = ensemble;
        cap_ens = cap_ens ? cap_ens * 2 : 1024;
        ensemble = alloue(cap_ens * n);
        nb_ens = 0;
        for (size_t i = 0; i < ancien; i++)
            if (vieux[i * n]) insere(vieux + i * n);
        free(vieux);
    }
    size_t i = hache(s, n) & (cap_ens - 1);
    while (ensemble[i * n]) {
        if (!memcmp(ensemble + i * n, s, n)) return 0;
        i = (i + 1) & (cap_ens - 1);
    }
    memcpy(ensemble + i * n, s, n);
    nb_ens++;
    return 1;
}

static void dessine(const char *sol);

static void solution(const int *lignes, int k) {
    int n = Z.nb_cases;
    char *s = alloue(n + 1), *img = alloue(n + 1), *canon = alloue(n + 1);
    for (int i = 0; i < k; i++) {
        const Placement *pl = &placements[lignes[i]];
        for (int c = 0; c < Z.pieces[pl->piece].n; c++) s[pl->cases[c]] = Z.pieces[pl->piece].nom;
    }
    /* forme canonique = la plus petite image par les symétries du plateau ;
     * on compte aussi combien d'images distinctes elle a (taille de l'orbite) */
    memcpy(canon, s, n);
    char images[48][128];
    int nb_images = 0;
    for (int t = 0; t < Z.nb_sym; t++) {
        for (int c = 0; c < n; c++) {
            char l = s[c];
            if (Z.sym_miroir[t]) l = l == Z.miroir_a ? Z.miroir_b : l == Z.miroir_b ? Z.miroir_a : l;
            img[Z.sym[t][c]] = l;
        }
        if (memcmp(img, canon, n) < 0) memcpy(canon, img, n);
        int nouveau = 1;
        for (int i = 0; i < nb_images && nouveau; i++) nouveau = memcmp(images[i], img, n) != 0;
        if (nouveau && n < 128) memcpy(images[nb_images++], img, n);
    }
    nb_brutes++;
    if (insere(canon)) {
        nb_total += nb_images;
        if (!premiere) { premiere = alloue(n + 1); memcpy(premiere, s, n); }
        if (a_afficher > 0) {
            a_afficher--;
            dessine(s);
            if (nb_gardees < 64) { gardees[nb_gardees] = alloue(n + 1); memcpy(gardees[nb_gardees++], s, n); }
        }
    }
    free(s); free(img); free(canon);
}

/* ================================================================ dessin */

static int couleur(char nom) {
    static const int pal[] = {196, 208, 226, 46, 51, 27, 129, 201, 94, 250, 22, 93, 166, 37};
    const char *ordre = "FILNPTUVWXYZAB";
    const char *q = strchr(ordre, nom);
    return pal[q ? q - ordre : 0];
}

static void dessine(const char *sol) {
    for (int y = 0; y < Z.H; y++) {
        for (int z = 0; z < Z.P; z++) {
            for (int x = 0; x < Z.L; x++) {
                int c = case_en(x, y, z);
                if (c < 0) fputs("  ", stdout);
                else if (couleurs) printf("\033[48;5;%dm\033[38;5;16m%c \033[0m", couleur(sol[c]), sol[c]);
                else printf("%c ", sol[c]);
            }
            if (z + 1 < Z.P) fputs("   ", stdout);
        }
        putchar('\n');
    }
    putchar('\n');
}


/* SVG : une case = un carré coloré, trait épais entre deux pièces différentes */
static const char *couleur_svg(char nom) {
    static const char *pal[] = {"#e8425a", "#f28c28", "#f2d02b", "#4cc76a", "#3cc7d6", "#3b6fe0",
                                "#8e44d8", "#e04fb0", "#9b6b43", "#c9c9c9", "#2e8b57", "#61008a", "#d2691e", "#20a0a0"};
    const char *ordre = "FILNPTUVWXYZAB";
    const char *q = strchr(ordre, nom);
    return pal[q ? q - ordre : 0];
}

static void ecrit_svg(const char *chemin) {
    FILE *f = fopen(chemin, "w");
    if (!f) meurt("impossible d'écrire le SVG");
    const int c = 24, marge = 12;
    int largeur_bloc = Z.P * Z.L * c + (Z.P - 1) * c, hauteur_bloc = Z.H * c;
    int par_ligne = nb_gardees < 4 ? nb_gardees : (Z.L > Z.H ? 2 : 4);
    if (par_ligne < 1) par_ligne = 1;
    int lignes = (nb_gardees + par_ligne - 1) / par_ligne;
    int W = par_ligne * (largeur_bloc + marge) + marge, H = lignes * (hauteur_bloc + marge) + marge;
    fprintf(f, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\" viewBox=\"0 0 %d %d\">\n", W, H, W, H);
    fprintf(f, "<rect width=\"%d\" height=\"%d\" fill=\"#0b0710\"/>\n", W, H);
    for (int k = 0; k < nb_gardees; k++) {
        const char *sol = gardees[k];
        int ox = marge + (k % par_ligne) * (largeur_bloc + marge), oy = marge + (k / par_ligne) * (hauteur_bloc + marge);
        for (int z = 0; z < Z.P; z++)
            for (int y = 0; y < Z.H; y++)
                for (int x = 0; x < Z.L; x++) {
                    int i = case_en(x, y, z);
                    if (i < 0) continue;
                    int px = ox + (z * (Z.L + 1) + x) * c, py = oy + y * c;
                    fprintf(f, "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" fill=\"%s\"/>\n", px, py, c, c, couleur_svg(sol[i]));
                    /* bords : là où la case voisine n'est pas de la même pièce */
                    int vx[4] = {1, -1, 0, 0}, vy[4] = {0, 0, 1, -1};
                    for (int d = 0; d < 4; d++) {
                        int j = case_en(x + vx[d], y + vy[d], z);
                        if (j >= 0 && sol[j] == sol[i]) continue;
                        int x1 = px + (vx[d] == 1 ? c : 0), y1 = py + (vy[d] == 1 ? c : 0);
                        int x2 = x1 + (vy[d] ? c : 0), y2 = y1 + (vx[d] ? c : 0);
                        fprintf(f, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#0b0710\" stroke-width=\"3\" stroke-linecap=\"square\"/>\n", x1, y1, x2, y2);
                    }
                }
    }
    fputs("</svg>\n", f);
    fclose(f);
}

/* ============================================================== puzzles */

static void pentaminos(void) {
    Z.plan = 1;
    ajoute_piece('F', ".##/##./.#.");
    ajoute_piece('I', "#####");
    ajoute_piece('L', "####/#...");
    ajoute_piece('N', "##../.###");
    ajoute_piece('P', "##/##/#.");
    ajoute_piece('T', "###/.#./.#.");
    ajoute_piece('U', "#.#/###");
    ajoute_piece('V', "#../#../###");
    ajoute_piece('W', "#../##./.##");
    ajoute_piece('X', ".#./###/.#.");
    ajoute_piece('Y', "####/.#..");
    ajoute_piece('Z', "##./.#./.##");
}

static void soma(void) {
    Z.plan = 0;
    ajoute_piece('V', "##/#.");
    ajoute_piece('L', "###/#..");
    ajoute_piece('T', "###/.#.");
    ajoute_piece('Z', ".##/##.");
    ajoute_piece('A', "##/.#|../.#"); /* vis à droite */
    ajoute_piece('B', "##/.#|#./.."); /* vis à gauche : image de A dans un miroir */
    ajoute_piece('P', "##/#.|#./..");  /* trépied */
    Z.miroir_a = 'A';
    Z.miroir_b = 'B';
}

static char *lit_fichier(const char *chemin) {
    FILE *f = fopen(chemin, "rb");
    if (!f) return NULL;
    size_t cap = 4096, n = 0;
    char *b = alloue(cap);
    int c;
    while ((c = fgetc(f)) != EOF) {
        if (n + 2 > cap) b = agrandit(b, cap *= 2);
        b[n++] = (char)c;
    }
    b[n] = 0;
    fclose(f);
    while (n && (b[n - 1] == '\n' || b[n - 1] == '\r')) b[--n] = 0;
    return b;
}

static char *rectangle(int l, int h) {
    char *t = alloue((size_t)(l + 1) * h + 1), *p = t;
    for (int y = 0; y < h; y++) { memset(p, '#', l); p += l; *p++ = '\n'; }
    return t;
}

static double maintenant(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main(int argc, char **argv) {
    const char *quoi = "6x10", *svg = NULL;
    int afficher = 1, briser = 1;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--afficher") && i + 1 < argc) afficher = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--sans-astuce")) briser = 0;
        else if (!strcmp(argv[i], "--svg") && i + 1 < argc) svg = argv[++i];
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--aide")) {
            puts("pavages [6x10|5x12|4x15|3x20|8x8|soma|fichier.txt] [--afficher N] [--svg F] [--sans-astuce]\n\n"
                 "  6x10 ... 3x20  les 12 pentaminos dans un rectangle\n"
                 "  8x8            un carré 8x8 avec un trou 2x2 au centre\n"
                 "  soma           le cube Soma (7 pièces, cube 3x3x3)\n"
                 "  fichier.txt    une forme de 60 cases : '#' à remplir, '.' ou espace = trou\n"
                 "  --afficher N   dessiner les N premières solutions distinctes (défaut 1)\n"
                 "  --svg F        les enregistrer aussi en image SVG\n"
                 "  --sans-astuce  ne pas restreindre la pièce X (plus lent, même résultat)");
            return 0;
        } else quoi = argv[i];
    }
    genere_mats();
    couleurs = isatty(1);
    char *texte;
    int prof = 1;
    int l, h;
    if (!strcmp(quoi, "soma")) {
        soma();
        texte = "###\n###\n###\n###\n###\n###\n###\n###\n###";
        prof = 3;
    } else {
        pentaminos();
        if (sscanf(quoi, "%dx%d", &h, &l) == 2 && strcmp(quoi, "8x8") && l * h == 60) texte = rectangle(l, h);
        else if (!strcmp(quoi, "8x8")) texte = "########\n########\n########\n###..###\n###..###\n########\n########\n########";
        else if (!(texte = lit_fichier(quoi))) meurt("puzzle inconnu ou fichier illisible (voir --aide)");
    }
    plateau_depuis_texte(texte, prof);
    int cubes = 0;
    for (int i = 0; i < Z.nb_pieces; i++) cubes += Z.pieces[i].n;
    if (cubes != Z.nb_cases) {
        fprintf(stderr, "pavages : la forme a %d cases mais les pièces en couvrent %d\n", Z.nb_cases, cubes);
        return 1;
    }
    int orient = 0;
    for (int i = 0; i < Z.nb_pieces; i++) { orientations(&Z.pieces[i]); orient += Z.pieces[i].nb_orient; }
    symetries();
    construit_matrice(briser);
    printf("%d pièces, %d orientations, %d cases, %d symétries du plateau -> matrice %d lignes × %d colonnes\n\n",
           Z.nb_pieces, orient, Z.nb_cases, Z.nb_sym, nb_placements, Z.nb_pieces + Z.nb_cases);
    a_afficher = afficher;
    double t0 = maintenant();
    cherche(0, solution);
    double dt = maintenant() - t0;
    if (!nb_ens) puts("aucune solution.");
    if (svg && nb_gardees) ecrit_svg(svg);
    printf("%zu solution(s) distincte(s) aux symétries près, %lld en tout\n", nb_ens, nb_total);
    printf("%.2f s, %.0f millions de mises à jour de liens (%.0f M/s)\n", dt, mises_a_jour / 1e6, mises_a_jour / 1e6 / (dt > 0 ? dt : 1));
    return 0;
}
