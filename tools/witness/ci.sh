#!/bin/sh
# ci.sh "<message>" files...
CAMR_HOME=${CAMR_HOME:-$(cd "$(dirname "$0")/../.." && pwd)}
cd $CAMR_HOME; m=$1; shift
git add "$@" && git commit -q -m "$m

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Claude-Session: ${CLAUDE_SESSION:-none}" 2>&1 | grep -v "unable to unlink"
$CAMR_HOME/tools/witness/sweep.sh > /dev/null; git log --oneline | head -1 | cut -c1-120; GIT_OPTIONAL_LOCKS=0 git status --short | grep -v '^??'
