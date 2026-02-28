#!/usr/bin/env bash
set -euo pipefail

# Atualiza o pipeline do VIPER em /opt/viper/src (ou no caminho definido em $REPO).
# Variaveis opcionais:
#   REPO   - caminho do repo (padrao: /opt/viper/src)
#   BRANCH - branch remota alvo (padrao: main)

REPO="${REPO:-/opt/viper/src}"
BRANCH="${BRANCH:-main}"

if [ ! -d "$REPO/.git" ]; then
  echo "Repositorio git nao encontrado em $REPO"
  exit 1
fi

cd "$REPO"

echo "Buscando atualizacoes em $REPO..."
git fetch --quiet origin

UPSTREAM="$(git rev-parse --abbrev-ref --symbolic-full-name @{u} 2>/dev/null || true)"
if [ -z "$UPSTREAM" ]; then
  UPSTREAM="origin/$BRANCH"
fi

LOCAL="$(git rev-parse HEAD)"
REMOTE="$(git rev-parse "$UPSTREAM")"

if [ "$LOCAL" = "$REMOTE" ]; then
  echo "Pipeline ja esta atualizado."
else
  echo "Atualizando a partir de $UPSTREAM..."
  git merge --ff-only "$UPSTREAM"
fi

# Espaco para comandos adicionais apos o update (preencha conforme necessario).
POST_HOOK="$REPO/update.d/post-update.sh"
if [ -f "$POST_HOOK" ]; then
  echo "Executando pos-atualizacao: $POST_HOOK"
  bash "$POST_HOOK"
fi

echo "Atualizacao concluida."
