#############
# relax functions
#############

from alphafold.relax import relax
from alphafold.common import protein

from simtk.openmm import app
from simtk.unit import nanometers, sqrt

# applied https://raw.githubusercontent.com/deepmind/alphafold/main/docker/openmm.patch
# to OpenMM 7.5.1 (see PR https://github.com/openmm/openmm/pull/3203)
# patch is licensed under CC-0
# OpenMM is licensed under MIT and LGPL
# fmt: off
def createDisulfideBonds(self, positions):
  def isCyx(res):
    names = [atom.name for atom in res._atoms]
    return 'SG' in names and 'HG' not in names
  # This function is used to prevent multiple di-sulfide bonds from being
  # assigned to a given atom.
  def isDisulfideBonded(atom):
    for b in self._bonds:
      if (atom in b and b[0].name == 'SG' and
        b[1].name == 'SG'):
        return True

    return False

  cyx = [res for res in self.residues() if res.name == 'CYS' and isCyx(res)]
  atomNames = [[atom.name for atom in res._atoms] for res in cyx]
  for i in range(len(cyx)):
    sg1 = cyx[i]._atoms[atomNames[i].index('SG')]
    pos1 = positions[sg1.index]
    candidate_distance, candidate_atom = 0.3*nanometers, None
    for j in range(i):
      sg2 = cyx[j]._atoms[atomNames[j].index('SG')]
      pos2 = positions[sg2.index]
      delta = [x-y for (x,y) in zip(pos1, pos2)]
      distance = sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2])
      if distance < candidate_distance and not isDisulfideBonded(sg2):
        candidate_distance = distance
        candidate_atom = sg2
    # Assign bond to closest pair.
    if candidate_atom:
      self.addBond(sg1, candidate_atom)
# fmt: on
app.Topology.createDisulfideBonds = createDisulfideBonds

def run_relax(pdb_filename=None, pdb_lines=None, pdb_obj=None, use_gpu=False):
  if pdb_obj is None:    
    if pdb_lines is None:
      pdb_lines = Path(pdb_filename).read_text()
    pdb_obj = protein.from_pdb_string(pdb_lines)
  
  amber_relaxer = relax.AmberRelaxation(
    max_iterations=0,
    tolerance=2.39,
    stiffness=10.0,
    exclude_residues=[],
    max_outer_iterations=3,
    use_gpu=use_gpu)
  
  relaxed_pdb_lines, _, _ = amber_relaxer.process(prot=pdb_obj)
  return relaxed_pdb_lines