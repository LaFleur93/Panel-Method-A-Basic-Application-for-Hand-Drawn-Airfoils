from tkinter import *
from tkinter import ttk
from tkinter import filedialog
import os
import matplotlib.pyplot as plt
from matplotlib import path
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import scipy.interpolate as interpolate
import numpy as np
from LineBuilder import *
from parametrization import parametric_airfoil
from panel_method_data import PanelMethod

class AirfoilPanelMethod:

	def __init__(self, window):

		#-------------------------------------
		# Main window
		#-------------------------------------
		self.window = window
		self.window.title('Panel Method, a Basic App for Hand-drawn Airfoils')
		self.window.geometry('425x615')
		self.window.resizable(False, False)

		#-------------------------------------
		# Left column - Airfoil geometry / Parameters
		#-------------------------------------

		lc_xpos = 10
		lc_width = .955

		#-------------------------------------
		# Message frame
		#-------------------------------------
		message_frame = LabelFrame(self.window, text = 'Message Box')
		message_frame.place(x = lc_xpos, y = 10, relwidth= lc_width, relheight = .1)

		self.message = Label(message_frame, text= 'No Message', fg = 'black', justify = "center", font = ('Helvetica', 10, 'bold'))
		self.message.place(x = lc_xpos + 10, y = 8)

		#-------------------------------------
		# Input parameters
		#-------------------------------------

		if_xpos = 145
		if_ypos = 2

		input_parameters_frame = LabelFrame(window, text = 'Input Parameters')
		input_parameters_frame.place(x = lc_xpos, y = 75, relwidth= .745, relheight = .2)

		Label(input_parameters_frame, text = 'Airfoil Chord (c)').place(x = lc_xpos, y = if_ypos)
		Label(input_parameters_frame, text = '[m]').place(x = lc_xpos + if_xpos, y = if_ypos)
		self.chord = Entry(input_parameters_frame, width= 7)
		self.chord.place(x = lc_xpos + if_xpos + 53, y = if_ypos)
		self.chord.insert(0,'1')
		Button(input_parameters_frame, text = 'SET', relief="solid", bg = '#CCCCCC', bd = 1, command = lambda: self.add_chord(self.chord)).place(x = lc_xpos + if_xpos + 105, y = if_ypos, relwidth= .14, relheight = .18)

		Label(input_parameters_frame, text = 'Airfoil Angle of Attack (α)').place(x = lc_xpos, y = if_ypos + 25*1)
		Label(input_parameters_frame, text = '[º]').place(x = lc_xpos + if_xpos, y = if_ypos + 25*1)
		self.AoA = Entry(input_parameters_frame, width= 7)
		self.AoA.place(x = lc_xpos + if_xpos + 53, y = if_ypos + 25*1)
		self.AoA.insert(0,'0')
		Button(input_parameters_frame, text = 'SET', relief="solid", bg = '#CCCCCC', bd = 1, command = lambda: self.add_AoA(self.AoA)).place(x = lc_xpos + if_xpos + 105, y = if_ypos + 25*1, relwidth= .14, relheight = .18)

		Label(input_parameters_frame, text = 'Free Stream Velocity (V∞)').place(x = lc_xpos, y = if_ypos + 25*2)
		Label(input_parameters_frame, text = '[m/s]').place(x = lc_xpos + if_xpos, y = if_ypos + 25*2)
		self.V_inf = Entry(input_parameters_frame, width= 7)
		self.V_inf.place(x = lc_xpos + if_xpos + 53, y = if_ypos + 25*2)
		self.V_inf.insert(0,'1')
		Button(input_parameters_frame, text = 'SET', relief="solid", bg = '#CCCCCC', bd = 1, command = lambda: self.add_velocity(self.V_inf)).place(x = lc_xpos + if_xpos + 105, y = if_ypos + 25*2, relwidth= .14, relheight = .18)
		
		Label(input_parameters_frame, text = 'Free Stream Density (ρ∞)').place(x = lc_xpos, y = if_ypos + 25*3)
		Label(input_parameters_frame, text = '[kg/m3]').place(x = lc_xpos + if_xpos, y = if_ypos + 25*3)
		self.rho_inf = Entry(input_parameters_frame, width= 7)
		self.rho_inf.place(x = lc_xpos + if_xpos + 53, y = if_ypos + 25*3)
		self.rho_inf.insert(0,'1.225')
		Button(input_parameters_frame, text = 'SET', relief="solid", bg = '#CCCCCC', bd = 1, command = lambda: self.add_density(self.rho_inf)).place(x = lc_xpos + if_xpos + 105, y = if_ypos + 25*3, relwidth= .14, relheight = .18)

		input_values_frame = LabelFrame(window, text = 'Values')
		input_values_frame.place(x = lc_xpos + if_xpos + 180, y = if_ypos + 73, relwidth= .19, relheight = .2)

		self.chord_box = Label(input_values_frame, text = 'c= None')
		self.chord_box.place(x = lc_xpos - 8, y = 0)
		self.AoA_box = Label(input_values_frame, text = 'α= None')
		self.AoA_box.place(x = lc_xpos - 8, y = 25*1)
		self.V_inf_box = Label(input_values_frame, text = 'V∞= None')
		self.V_inf_box.place(x = lc_xpos - 8, y = 25*2)
		self.rho_inf_box = Label(input_values_frame, text = 'ρ∞= None')
		self.rho_inf_box.place(x = lc_xpos - 8, y = 25*3)

		#-------------------------------------
		# Main Buttons (Draw / Load / Calculate)
		#-------------------------------------

		img_ypos = 205
		button_font = ('Raleway', 10, 'bold')

		self.img = ones((100,100,3))
		draw_button = Button(window, text = 'Draw', font = button_font, bg = '#007F7F', fg = 'white', command = lambda: self.draw_airfoil(self.img))
		draw_button.place(x = lc_xpos, y = img_ypos, relwidth= 0.23, relheight = .05)

		load_button = Button(window, text = 'Load', font = button_font, bg = '#007F7F', fg = 'white', command = lambda: self.load_airfoil())
		load_button.place(x = lc_xpos + 102, y = img_ypos, relwidth= 0.23, relheight = .05)
		
		self.nodes = []
		self.AoA_value = ''
		self.V_inf_value = ''
		self.rho_inf_value = ''
		self.chord_value = ''

		self.calculate_button = Button(window, text = 'Click To Calculate', relief="solid", bg = '#CCCCCC', bd = 1, command = lambda: self.calculate(self.nodes, self.AoA_value, self.V_inf_value, self.rho_inf_value, self.chord_value))
		self.calculate_button.place(x = lc_xpos + 205 , y = img_ypos, relwidth= 0.47, relheight = .05)

		#-------------------------------------
		# Drawn airfoil
		#-------------------------------------

		self.parametrized_airfoil_graphic = Frame(window, bg='#C0C0C0', bd=1.5)
		self.parametrized_airfoil_graphic.place(x = lc_xpos, y = img_ypos + 40, relwidth = lc_width)
		self.parametrized_figure = plt.Figure(figsize=(0.05,2.5), dpi=80)
		self.parametrized_axes = self.parametrized_figure.add_subplot(111)
		self.parametrized_axes.set_ylabel("x/c vs. y/c")
		self.parametrized_axes.grid(True)
		self.parametrized_line = FigureCanvasTkAgg(self.parametrized_figure, self.parametrized_airfoil_graphic)
		self.parametrized_line.get_tk_widget().pack(side=LEFT, fill=BOTH, expand=1)
        
		#-------------------------------------
		# Output Data
		#-------------------------------------
	
		output_data_frame = LabelFrame(window, bg = 'white')
		output_data_frame.place(x = lc_xpos, y = 455, relwidth= lc_width, relheight = .05)
		output_data_message = Label(output_data_frame, text = 'Output Data', fg = 'black', bg = 'white', justify = "center", font = ('Helvetica', 10, 'bold'))
		output_data_message.place(x = 5, y = 2)

		#-------------------------------------
		# Output parameters (Cl, Lift, etc.)
		#-------------------------------------

		output_xpos = 177 # Second Column
		output_ypos = 497
		output_font = ('Raleway', 10, 'normal')

		Label(window, text = 'Lift Coefficient (Cl)', font = output_font).place(x = lc_xpos, y = output_ypos)
		self.lift_coefficient = Label(window, text = 'Cl= None', font = output_font)
		self.lift_coefficient.place(x = output_xpos, y = output_ypos)

		Label(window, text = 'Airfoil Lift (L)', font = output_font).place(x = lc_xpos, y = output_ypos + 25*1)
		self.airfoil_lift = Label(window, text = 'L= None', font = output_font)
		self.airfoil_lift.place(x = output_xpos, y = output_ypos + 25*1)

		Label(window, text = 'Moment Coefficient (Cm)', font = output_font).place(x = lc_xpos, y = output_ypos + 25*2)
		self.moment_coefficient = Label(window, text = 'Cm= None', font = output_font)
		self.moment_coefficient.place(x = output_xpos, y = output_ypos + 25*2)

		Label(window, text = 'Airfoil Moment (M)', font = output_font).place(x = lc_xpos, y = output_ypos + 25*3)
		self.airfoil_moment = Label(window, text = 'M= None', font = output_font)
		self.airfoil_moment.place(x = output_xpos, y = output_ypos + 25*3)

		#-------------------------------------
		# Graphics buttons
		#-------------------------------------

		streamplot_button = Button(window, text = 'Stream Plot', font = button_font, bg = '#007F7F', fg = 'white', command = self.streamplot_graphic)
		streamplot_button.place(x = 275, y = 493, relwidth = 0.33)

		pressure_coefficient_button = Button(window, text = 'Cp Over Airfoil', font = button_font, bg = '#007F7F', fg = 'white', command = self.cp_airfoil_graphic)
		pressure_coefficient_button.place(x = 275, y = 493 + 28*1, relwidth = 0.33)

		pressure_contour_button = Button(window, text = 'Cp Contour Plot', font = button_font, bg = '#007F7F', fg = 'white', command = self.cp_contour_graphic)
		pressure_contour_button.place(x = 275, y = 493 + 28*2, relwidth = 0.33)

		about_button = Button(window, text = 'About', font = button_font, bg = '#007F7F', fg = 'white', command = self.about)
		about_button.place(x = 275, y = 493 + 28*3, relwidth = 0.33)

	def add_chord(self, parameter):
		if self.validate_input(parameter):
			self.chord_box['text'] = f'c= {parameter.get()} m'
			self.chord_value = float(parameter.get())
	
	def add_AoA(self, parameter):
		if self.validate_input(parameter, True):
			self.AoA_box['text'] = f'α= {parameter.get()}º'
			self.AoA_value = float(parameter.get())

	def add_velocity(self, parameter):
		if self.validate_input(parameter):
			self.V_inf_box['text'] = f'V∞= {parameter.get()} m/s'
			self.V_inf_value = float(parameter.get())

	def add_density(self, parameter):
		if self.validate_input(parameter):
			self.rho_inf_box['text'] = f'ρ∞= {parameter.get()}'
			self.rho_inf_value = float(parameter.get())

	def validate_input(self, parameter, AoA = False):
		try:
			p = float(parameter.get())
			if p > 0 and AoA == False:
				self.message['text'] = 'Parameter succesfully updated'
				self.message['fg'] = 'green'
				return True
			elif p >= 0 and AoA == True:
				self.message['text'] = 'Parameter succesfully updated'
				self.message['fg'] = 'green'
				return True
			else:
				self.message['text'] = 'Parameter should be a positive number'
				self.message['fg'] = 'red'
				return False
		except:
			self.message['text'] = 'Parameter should be a number'
			self.message['fg'] = 'red'
			return False

	def draw_airfoil(self, data):
		x_min, x_max = -.1, 1.1
		y_min, y_max = -.25, .25
		fig = plt.figure()
		ax = fig.add_subplot(111)
		#ax.set_xticks(np.arange(-0.2, 1.2, 0.1))
		ax.set_yticks(np.arange(-0.5, 0.5, 0.05))
		ax.set_axisbelow(True)
		ax.grid(True, ls = '--')

		currentAxis = plt.gca()
		currentAxis.add_patch(Rectangle((0,-0.1), 1, 0.2, alpha=1, ec = 'gray', fc = 'None'))

		line = ax.imshow(data)
		ax.set_xlim(x_min, x_max)
		ax.set_ylim(y_min, y_max)
		plt.title("Start from Trailing Edge and continue in Clockwise Direction")
		plt.xlabel("Distance - x/c")
		plt.ylabel("Distance - y/c")

		draw_airfoil = LineBuilder(line,ax)
		plt.show()
		
		if draw_airfoil.nodes != []:
			self.airfoil = draw_airfoil
			self.parametric_data = parametric_airfoil(self.airfoil.nodes)
			self.parametrized_axes.clear()
			self.parametrized_axes.plot(self.parametric_data[0][0], self.parametric_data[0][1], linewidth=2, color = 'gray')
			self.parametrized_axes.plot(self.parametric_data[1][:,0], self.parametric_data[1][:,1], 'ob', color = 'gray', mec = 'black')
			self.parametrized_axes.set_xlim(x_min, x_max)
			self.parametrized_axes.set_ylim(y_min, y_max)
			self.parametrized_axes.set_ylabel("x/c vs. y/c")
			self.parametrized_axes.grid(True)
			self.parametrized_line.draw()
			
			self.nodes = []
			for i in range(len(self.parametric_data[0][0])):
				self.nodes.append((self.parametric_data[0][0][i], self.parametric_data[0][1][i]))

	def load_airfoil(self):
		root = os.path.join(os.getcwd(), 'SampleAirfoils')
		file = filedialog.askopenfilename(initialdir = root)
		if file:
			airfoil_data = np.loadtxt(file, delimiter = '   ', skiprows = 1)
			nodes = np.flipud(airfoil_data)
			nodes = nodes[:-1]
			first_node = nodes[0]
			nodes = np.vstack([nodes, first_node])
			x_min, x_max = -.1, 1.1
			y_min, y_max = -.25, .25
			self.parametrized_axes.clear()
			self.parametrized_axes.plot(airfoil_data[:,0], airfoil_data[:,1], linewidth=2, color = 'gray')
			self.parametrized_axes.plot(airfoil_data[:,0], airfoil_data[:,1], 'ob', color = 'gray', mec = 'black')
			self.parametrized_axes.set_xlim(x_min, x_max)
			self.parametrized_axes.set_ylim(y_min, y_max)
			self.parametrized_axes.set_ylabel("x/c vs. y/c")
			self.parametrized_axes.grid(True)
			self.parametrized_line.draw()
			
			self.nodes = []

			for i in range(len(nodes)):
				self.nodes.append((nodes[:,0][i], nodes[:,1][i]))

	def calculate(self, nodes, angle, velocity, density, chord):
		if self.validate_data():
			
			self.data = PanelMethod(nodes, self.AoA_value, self.V_inf_value, self.rho_inf_value, self.chord_value)
			self.xps = self.data[0]
			self.yps = self.data[1]
			self.U = self.data[2]
			self.V = self.data[3]
			self.Cp_XY = self.data[4]
			self.x = self.data[5]
			self.Cp = self.data[6]
			self.Cl = round(self.data[7], 3)
			self.L = round(self.data[8], 2)
			self.Cm = round(self.data[9], 3)
			self.M = round(self.data[10], 2)
			
			self.lift_coefficient['text'] = f'Cl= {self.Cl}'
			self.airfoil_lift['text'] = f'L= {self.L} N'
			self.moment_coefficient['text'] = f'Cm= {self.Cm}'
			self.airfoil_moment['text'] = f'M= {self.M} Nm'

			self.message['text'] = 'Data succesfully calculated'
			self.message['fg'] = 'green'
		
	def validate_data(self):
		try:
			float(self.chord_value)
			float(self.AoA_value)
			float(self.V_inf_value)
			float(self.rho_inf_value)
			if self.nodes != []:
				return True
			else:
				self.message['text'] = 'Please, Draw or Load Airfoil'
				self.message['fg'] = 'red'
				return False
		except:
			self.message['text'] = 'Please, Check Inputs'
			self.message['fg'] = 'red'
			return False

	def streamplot_graphic(self):
		if self.validate_calculated_data():
			fig, strm_ax = plt.subplots()
			strm = strm_ax.streamplot(self.xps, self.yps, self.U, self.V, color=(self.U**2+self.V**2)**0.5, linewidth=2, cmap=plt.cm.autumn)
			fig.colorbar(strm.lines)
			strm_ax.set_title('Velocity Over Airfoil [m/s]')
			strm_ax.set_xlabel('Distance - x/c')
			strm_ax.set_ylabel('Distance - y/c')
			strm_ax.set_axisbelow(True)
			x = [node[0] for node in self.nodes]
			y = [node[1] for node in self.nodes]
			plt.fill_between(x, y, color = 'gray')
			plt.show()

	def cp_airfoil_graphic(self):
		if self.validate_calculated_data():
			plt.scatter(self.x, self.Cp)
			plt.plot(self.x, self.Cp)
			ax = plt.gca()
			ax.set_title('Pressure Coefficient over Airfoil')
			ax.set_xlabel('Distance - x/c')
			ax.set_ylabel('Pressure Coefficient - Cp')
			ax.grid(True)
			ax.set_ylim(ax.get_ylim()[::-1])
			plt.show()

	def cp_contour_graphic(self):
		if self.validate_calculated_data():
			fig = plt.figure()
			plt.contourf(self.xps, self.yps, self.Cp_XY,500,cmap='jet')
			plt.title('Pressure Coefficient Contour Plot')
			plt.xlabel('Distance - x/c')
			plt.ylabel('Pressure Coefficient - Cp')
			x = [node[0] for node in self.nodes]
			y = [node[1] for node in self.nodes]
			plt.fill_between(x, y, color = 'gray')
			plt.colorbar()
			plt.show()

	def about(self):
		path = 'About.pdf'
		os.system(path)

	def validate_calculated_data(self):
		try:
			self.data
			return True
		except:
			self.message['text'] = 'No Data Calculated. Plase Check Inputs and/or Airfoil'
			self.message['fg'] = 'red'
			return False

if __name__ == "__main__":
	window = Tk()
	application = AirfoilPanelMethod(window)
	window.mainloop()