import typing as T


def _parse_indices(val: T.Optional[str]) -> T.Optional[list]:
    if not val:
        return None
    parts = [p.strip() for p in val.split(",") if p.strip()]
    ints = []
    for p in parts:
        try:
            ints.append(int(p))
        except ValueError:
            continue
    return ints if ints else None


def _map_indices_to_aligned(
    aligned_seq: str, indices: T.Optional[T.List[int]]
) -> T.List[int]:
    if not indices:
        return []
    mapped: T.List[int] = []
    ungapped_counter = 0
    pos_map = {}
    for aligned_i, ch in enumerate(aligned_seq, start=1):
        if ch != "-":
            ungapped_counter += 1
            pos_map[ungapped_counter] = aligned_i
    for idx in indices:
        if 1 <= idx <= ungapped_counter:
            mapped.append(pos_map[idx])
    return mapped
