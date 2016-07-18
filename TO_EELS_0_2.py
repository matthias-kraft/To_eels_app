import Tkinter as tk
from TO_EELS_layout_0_2 import TO_layout
from App_position_variables_0_2 import *

# GUI main frame
root = tk.Tk()
root.minsize(width,height)
root.configure(background='black')
root.wm_title('EELS with Transformation optics')

TO_app = TO_layout(root)
root.mainloop()


