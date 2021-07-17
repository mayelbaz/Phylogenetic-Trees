from tkinter import *
from tkinter.ttk import Progressbar

import main
from PIL import ImageTk, Image
from threading import Thread


class GUI:
    # initialzing all the main components
    def __init__(self, master):
        self.progress_bar_appear = False
        self.show_alignment_btn = False
        self.show_tree_btn = False
        self.master = master
        master.title("Phylogenetic tree creator")

        self.control_frame = Frame(master, relief=RAISED, borderwidth=1, bg ="PaleVioletRed1")
        self.control_frame.pack(fill=BOTH, expand=True, side = TOP)

        self.log_frame = Frame(master, relief=RAISED, borderwidth=1, bg ="PaleVioletRed1")
        self.log_frame.pack(fill=BOTH, expand=True, side = LEFT)

        self.stats_frame = Frame(master, relief=RAISED, borderwidth=1, bg ="misty rose")
        self.stats_frame.pack(fill=BOTH, expand=True, side = RIGHT)

        # adding the buttons
        #control frame
        self.start_button = Button(self.control_frame, text="      GO!         ", command=self.start)
        self.start_button.pack(side=RIGHT,padx=10)

        self.close_button = Button(self.control_frame, text="  Close  ", command=self.quit)
        self.close_button.pack(side=LEFT)

        # adding the canvas
        #log frame
        self.canvas = Canvas(self.log_frame, width=600, height=600)
        self.canvas.pack()

        # adding the genes and side bar lists
        self.empty = Text(self.stats_frame, bg="MISTY ROSE3", width = 40, height = 2)
        self.empty.pack()

        self.gene_choice_option = StringVar(master)
        self.gene_choice_option.set("Choose a Gene")  # default value

        self.gene_choice = OptionMenu(self.stats_frame, self.gene_choice_option, "P53", "APOE", "BRCA1", "VEGFA", command=self.show_animals)
        self.gene_choice.pack()

        self.animal_choice = Listbox(self.stats_frame, selectmode="multiple")
        self.animal_choice.pack(expand=NO, fill="both")

        self.animal_text = Text(self.stats_frame, bg="MISTY ROSE3", width=40, height=2)
        #self.animal_text.configure(state="disabled")
        self.animal_text.pack()

        # add a schedualer every 10 MS to check if should show tree
        root.after(10,self.show_tree)

    # shows the animals after choosing the gene
    def show_animals(self,gene):
        self.animal_text.delete('1.0', END)
        self.animal_text.insert("1.0", "Choose 3-10 organisms")

        self.animal_choice.delete(0, 'end')

        x = main.get_animals_for_gene(gene)
        i = 0
        for each_item in x:
            self.animal_choice.insert(END, each_item)

            # coloring alternative lines of listbox
            self.animal_choice.itemconfig(i,
                            bg="THISTLE1" if i % 2 == 0 else "HotPink2")
            i=i+1
    # starts the process
    def start(self):
        gene = self.gene_choice_option.get()
        animals = [self.animal_choice.get(idx) for idx in self.animal_choice.curselection()]
        if len(animals) < 3:
            self.empty.insert("1.0", "ERROR: Please choose 3 to 10 animals")
            self.empty.configure(state="disabled")
            return

        self.empty.configure(state="normal")
        self.empty.delete('1.0', END)
        self.empty.configure(state="disabled")

        self.thread = Thread(target=main.start_program, args=(gene,animals))
        self.thread.start()
        if not self.progress_bar_appear:
            self.progress_bar = Progressbar(self.stats_frame,orient=HORIZONTAL,length=250,  mode='determinate')
            #self.animal_text.configure(state="disabled")
            self.progress_bar.pack()
            self.progress_bar_appear = True
        root.after(1,self.bar)
        root.after(10,self.show_tree)

    # updates the status bar
    def bar(self):
        self.progress_bar['value'] = main.job_done_precentage
        root.after(1,self.bar)

    # stops the program
    def quit(self):
        self.master.quit()

    # shows the alignment after clicking the button
    def show_alignment(self):
        if main.alignment_file_name != "":
            with open(main.alignment_file_name,"r") as f:
                self.log_frame.pack_forget()
                self.log_frame = Frame(self.master, relief=RAISED, borderwidth=1, bg="PaleVioletRed1")
                self.log_frame.pack(fill=BOTH, expand=True, side=LEFT)

                text = f.read()
                self.alignment_txt = Text(self.log_frame, bg="MISTY ROSE3", width=80, height=40)
                # self.alignment_txt.configure(state="disabled")
                self.alignment_txt.insert('1.0',text)
                self.alignment_txt.pack()
    # shows the finished tree
    def show_tree_final(self):
        if main.photo_ready != "":
            self.log_frame.pack_forget()
            self.log_frame = Frame(self.master, relief=RAISED, borderwidth=1, bg="PaleVioletRed1")
            self.log_frame.pack(fill=BOTH, expand=True, side=LEFT)
            self.canvas = Canvas(self.log_frame, width=600, height=600)
            self.canvas.pack()
            self.img = ImageTk.PhotoImage(Image.open(main.saved_file_name))
            self.canvas.create_image(20, 20, anchor=NW, image=self.img)

    # schedualed function to show the tree when ready
    def show_tree(self):
        if main.photo_ready:
            self.show_tree_final()
            if not self.show_alignment_btn:
                self.show_alignment = Button(self.control_frame, text="Show Alignment", command=self.show_alignment)
                self.show_alignment.pack(side=TOP, padx=10)
                self.show_alignment_btn = True
            if not self.show_tree_btn:
                self.show_tree_actual_btn = Button(self.control_frame, text="Show Tree", command=self.show_tree_final)
                self.show_tree_actual_btn.pack(side=TOP, padx=10)
                self.show_tree_btn = True
            main.photo_ready = False
        else:
            root.after(10,self.show_tree)

if __name__ == '__main__':
    #starting the GUI
    root = Tk()
    my_gui = GUI(root)
    root.mainloop()
