"""Combinatorial free-variable planning, unchanged from v1.11."""


def plan(vw):
    known = set()
    order = []
    free = []
    pending = set(vw)

    def ready(v):
        word, es = vw[v]
        d = len(word)
        ks = [k for k in range(d) if es[k] in known]
        if d == 3:
            return []
        if d == 4 and ks:
            return ks[:1]
        if d == 5:
            for k in range(5):
                if es[k] in known and es[(k + 1) % 5] in known:
                    return [k, (k + 1) % 5]
        return None

    while pending:
        progressed = False
        for v in sorted(pending, key=lambda v: len(vw[v][0])):
            ks = ready(v)
            if ks is not None:
                order.append((v, ks))
                known |= set(vw[v][1])
                pending.discard(v)
                progressed = True
        if not progressed:
            cand = []
            for v in pending:
                word, es = vw[v]
                d = len(word)
                for k in range(d):
                    if es[k] in known:
                        continue
                    gain = (
                        1
                        if d == 4
                        else (
                            1
                            if (es[(k + 1) % d] in known or es[(k - 1) % d] in known)
                            else 0
                        )
                    )
                    cand.append((-gain, len(word), k, v))
            cand.sort()
            _, _, k, v = cand[0]
            free.append(vw[v][1][k])
            known.add(vw[v][1][k])
    return order, free
