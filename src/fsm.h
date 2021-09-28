#ifndef FSM_H
#define FSM_H

#include "state.h"

struct __gff_target;
struct gff_aux;
struct gff_tok;

static inline void fsm_init(enum state *state) { *state = STATE_BEGIN; }

enum state fsm_next(enum state state, struct gff_tok *tok, struct gff_aux *aux,
                    struct __gff_target *tgt);

char const *fsm_name(enum state state);

#endif
