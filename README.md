 # Processing and analysis of the sound in Julia

## Authors:
1. Jaworek Klaudia
2. Jelito Natalia
3. Spik Urszula
4. Szymkowiak Magdalena

## Technologies:
**Julia** - programming language, version 1.5.3 is recommended.

**OffsetArrays** - package providing arrays.

**WAV** - package that enable working with sound files.

**FFTW** - package used in processing sounds.

**Plots** - package for visualisation.

**Gtk** - package providing Graphical User Interface.

 ## Short description

 This project was created during the Mathematical Packages course by students from Faculty of Pure and Applied Mathematics at the Wroclaw University of Science and Technolody. 
 It is used for processing and analysis of the sound.Grapgical User Interface gives all of the funcionalities written in this project.

 ### Here are the functionalities:
 * changing frequency,
 * denoising based on the DSP,
 * removing frequencies,
 * lowpass and highpass,
 * changing speed,
 * changing volume,
 * cutting time,
 * smoothig signal.

## Some pictures showing our GUI:
![](/.png)



## Ho to run this program?
1. Clone the project to your directory: `git clone https://github.com/Klaudia226/Fourier-project.git`.
2. In Julia's terminal type:
   `using Pkg; Pkg.add(["OffsetArrays","WAV", "FFTW", "Plots","Gtk"])`
3. Open `appGUI.jl` in Visual Studio Code and click *Julia: Execude file in REPL* to run.