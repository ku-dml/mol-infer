# some convenient scripts
from collections import deque
import subprocess, time

def _read_clean_lines(path):
	q = []
	with open(path, 'r') as f:
		for raw in f:
			s = raw.strip()
			if not s or s.startswith("#"):
				continue
			q.append(s)
	return deque(q)

def _pop_line(Q):
	try:
		return Q.popleft()
	except IndexError:
		raise ValueError("Unexpected EOF while parsing seed-graph file.")

def _pop_int(Q):
	s = _pop_line(Q)
	try:
		return int(s)
	except ValueError:
		raise ValueError(f"Excepted integer line, got: {s!r}.")

def _pop_float(Q):
	s = _pop_line(Q)
	try:
		return float(s)
	except ValueError:
		raise ValueError(f"Excepted integer line, got: {s!r}.")

def _pop_ints(Q):
	s = _pop_line(Q)
	try:
		return list(map(int, s.split()))
	except ValueError:
		raise ValueError(f"Expected space-separated ints, got: {s!r}")

def _pop_floats(Q):
	s = _pop_line(Q)
	try:
		return list(map(float, s.split()))
	except ValueError:
		raise ValueError(f"Expected space-separated ints, got: {s!r}")

def _pop_strs(Q):
	s = _pop_line(Q)
	try:
		return list(s.split())
	except ValueError:
		raise ValueError(f"Expected space-separated strings, got: {s!r}")

def run_with_retry(cmd, retries=10, delay=2):
	for attempt in range(1 + retries + 1):
		try:
			subprocess.run(
				cmd,
				stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
				check=True, text=True
			)
			return
		except subprocess.CalledProcessError as e:
			print(f"[WARN] fv_gen failed (attempt {attempt}/{retries}): {e}")
			if attempt < retries:
				time.sleep(delay) 
			else:
				raise  

def _map_optional_list(x, mapping):
    if isinstance(x, (list, tuple)):
        return [mapping[v] for v in x]
    else:
        return mapping[x]

