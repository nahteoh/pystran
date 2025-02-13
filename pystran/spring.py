"""
Define spring mechanical quantities.
"""

from numpy import outer
from pystran import assemble
from pystran import freedoms


def assemble_stiffness(Kg, j):
    """
    Assemble the stiffness matrix of all the springs at joint `j`.

    - `Kg` is the global stiffness matrix,
    - `j` = the joint (obtained from the model as `m["joints"][jid]`).
    """
    dim = len(j["coordinates"])
    if "extension" in j["springs"]:
        for sd in j["springs"]["extension"].values():
            dof = [j["dof"][k] for k in freedoms.translation_dofs(dim)]
            coefficient, direction = sd["coefficient"], sd["direction"]
            k = coefficient * outer(direction, direction)
            Kg = assemble.assemble(Kg, dof, k)
    if "moment" in j["springs"]:
        for sd in j["springs"]["moment"].values():
            dof = [j["dof"][k] for k in freedoms.rotation_dofs(dim)]
            coefficient, direction = sd["coefficient"], sd["direction"]
            k = coefficient * outer(direction, direction)
            Kg = assemble.assemble(Kg, dof, k)
