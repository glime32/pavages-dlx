#!/usr/bin/env bash
# Compare les décomptes aux valeurs connues (quelques secondes en tout).
cd "$(dirname "$0")"
make -s pavages || exit 1
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
ok=0 total=0
verifie() { # description, attendu (morceau de la sortie), arguments de pavages...
  local desc="$1" attendu="$2"; shift 2
  total=$((total + 1))
  if timeout 20 ./pavages "$@" 2>&1 | grep -qF "$attendu"; then ok=$((ok + 1)); echo "  ok    $desc"
  else echo "  ÉCHEC $desc (attendu : $attendu)"; fi
}
verifie "3×20 : 2 pavages, 8 avec les symétries"      "2 solution(s) distincte(s) aux symétries près, 8 en tout" 3x20
verifie "4×15 : 368 pavages"                          "368 solution(s) distincte(s) aux symétries près, 1472 en tout" 4x15
verifie "8×8 troué au centre : 65 pavages"            "65 solution(s) distincte(s) aux symétries près, 520 en tout" 8x8
verifie "cube Soma : 240 solutions"                   "240 solution(s) distincte(s) aux symétries près, 11520 en tout" soma
verifie "cœur : 55 pavages"                           "55 solution(s) distincte(s) aux symétries près, 110 en tout" formes/coeur.txt
verifie "sans l'astuce sur X : même résultat"         "2 solution(s) distincte(s) aux symétries près, 8 en tout" 3x20 --sans-astuce
printf '###\n' > "$tmp/petit.txt"
verifie "forme qui n'a pas 60 cases : refusée"        "la forme a 3 cases" "$tmp/petit.txt"

total=$((total + 1))
if ./pavages 3x20 --afficher 2 --svg "$tmp/t.svg" > /dev/null && grep -q '^<svg ' "$tmp/t.svg" && tail -1 "$tmp/t.svg" | grep -q '</svg>'; then
  ok=$((ok + 1)); echo "  ok    --svg écrit un SVG complet"
else echo "  ÉCHEC --svg écrit un SVG complet"; fi

echo "$ok / $total"
[ "$ok" = "$total" ]
