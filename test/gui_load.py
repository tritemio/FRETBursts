from guidata.dataset.datatypes import DataSet
from guidata.dataset.dataitems import (FloatItem, IntItem, BoolItem, ChoiceItem,
                             MultipleChoiceItem, ImageChoiceItem, FilesOpenItem,
                             StringItem, TextItem, ColorItem, FileSaveItem,
                             FileOpenItem, DirectoryItem, FloatArrayItem,
                             DateItem, DateTimeItem)

SHOW=True

class GuiLoadData(DataSet):
    """
    Load a measurement file.
    """
    dir = DirectoryItem("Directory", "/home/anto/")
    fname = FileOpenItem("Open file", ("dat",), "")
    swap_DA = BoolItem("swap DA", "Swap Donor and Acceptor channels")

    
if __name__ == "__main__":
    # Create QApplication
    import guidata
    _app = guidata.qapplication()

    e = GuiLoadData()

