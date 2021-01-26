
from astroboi_bio_tools.ToolLogicPrep import ToolLogicPreps
class LogicPreps(ToolLogicPreps):
    def make_cell_id(self, brcd_arr, deli_str):
        result_list = []
        for brcd1 in brcd_arr:
            for brcd2 in brcd_arr:
                result_list.append(brcd1 + deli_str + brcd2)

        return result_list