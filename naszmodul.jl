"""
Informacje:
aby użyć modułu piszemy include("naszmodul.jl") (tak jak te wszystkie usingi)
A później używamy tych funkcji robiąc np. NaszModul.dftfreq(coś tam coś tam)
Nazwa NaszModul jest oczywiście dość robocza... xd
"""

"""Moduł z funkcjami do przetwarzania dźwięku"""
module NaszModul

using OffsetArrays #nie wiem czy wszystkie usingi mam
using WAV
using FFTW

"""Funkcja zwacająca macierz DFT"""
function dft(signal)
    N = length(signal)
    zeta_powers = OffsetArray([exp( -2π *  im * n / N) for n in 0:(N-1) ], 0:(N-1))
    [  sum( signal[n + 1] * zeta_powers[(n * f) % N] for n in 0:(N-1)   ) for f in 0:(N-1) ]
end

"""Funkcja zwracająca częstotliwość próbkowania DFT"""
function dftfreq(n, fs)
    if n >= 1
        f = zeros(n)
        if iseven(n)
            positive = collect(0:((n/2)-1))
            negative = collect((-n/2):1:(-1))
        else isodd(N)
            positive = collect(0:((n-1)/2))
            negative = collect((-(n-1)/2):1:(-1))
        end
        len = length(positive)
        f[1:len] = positive
        f[(len + 1):end] = negative
        return f*fs/n
    else
        throw(ArgumentError("n powinno być większe równe 1"))
    end
end

"""Odwrotna funkcja DFT"""
function idft(signal)
    N = length(signal)
    zeta_powers = OffsetArray([exp( 2π *  im * n / N) for n in 0:(N-1) ], 0:(N-1))
    [  (1/N)*sum( signal[n + 1] * zeta_powers[(n * f) % N] for n in 0:(N-1)   ) for f in 0:(N-1) ]
end

"""Funkcja wycinająca mniejsze wielkości(tzn. cichsze, na podstawie PSD)"""
function denoising_high(sound, fs, N,cut_point)
    tstep=1/fs
    time=LinRange(0, (N-1)/fs, N)
    freq_step = fs/N
    freq = LinRange(0, (N-1)*freq_step, N)
    
    fhat=fft(sound)
    PSD = real.(fhat .* conj(fhat)/N) #power spectral density
    
    indices=[i > cut_point for i in PSD]
    PSD_clean= indices .* PSD
    fhat_clean = indices .* fhat
    
    return time, PSD_clean, fhat_clean, freq
end

"""Funkcja wycinająca większe wielkości(tzn. głośniejsze, na podstawie PSD)"""
function denoising_low(sound, fs, N,cut_point)
    tstep=1/fs
    time=LinRange(0, (N-1)/fs, N) 
    freq_step = fs/N
    freq = LinRange(0, (N-1)/fs, N)
    
    fhat=fft(sound)
    PSD = real.(fhat .* conj(fhat)/N)
    
    indices=[i < cut_point for i in PSD]
    PSD_clean= indices .* PSD
    fhat_clean = indices .* fhat

    return time, PSD_clean, fhat_clean, freq
end

"""Funkcja wycinająca określoną częstotliwość"""
function remove_frequency(sound, freq_to_remove, dt = 0.001, start = 0, stop = 220.5)
    t = LinRange(start, stop, Int((stop - start)/dt))
    n = length(t) 
    fhat = real.(fft(sound[1][:, 1])) 
    PSD = fhat .* conj(fhat)/n 
    freq = 1/(dt * n) .* [i for i in 1:n] 
    
    indices = [i*dt*n != freq_to_remove for i in freq] 
    PSD_clean = indices .* PSD 
    fhat_clean = indices .* fhat 
    real(ifft(fhat_clean)), freq, PSD_clean, fhat_clean
end

"""Funkcja do przyspieszania/zwalniania tempa odtwarzania"""
function change_speed(sound, fs, speed, file_name)
    wavwrite(sound, file_name, Fs = fs * speed)
    sound_speeded, fs_speeded = wavread(file_name)
    return sound_speeded, fs_speeded
end

"""Funkcja do zmiany głośniści (ściszanie i zgłaśnianie)"""
function change_volume(sound, volume)
    sound_fft = fft(sound)
    new_sound_fft = volume .* sound_fft
    new_sound = ifft(new_sound_fft)
    return new_sound
end

"""Funkcja do przycinania czasu nagrania"""
function cutting_time(sound,fs,start=0,stop=((length(sound[:,1])-1)/fs))
    N=length(sound[:,1]) 
    time=LinRange(0,(N-1)/fs,N)
    
    if start==0
        new_sound=sound[:,1][1:floor(Int,stop/time[end]*N)]
    else 
        new_sound=sound[:,1][floor(Int,start/time[end]*N):floor(Int,stop/time[end]*N)]
    end
        
    N_new=length(new_sound)
    time_new=LinRange(0,(N_new-1)/fs,N_new)
    
    return time_new, new_sound
end

"""Oszacowanie wartości w zadanym punkcie za pomocą regresji lokalnej"""
function loess(index, X,Y,m)
    if index <= m
        xs = X[1:(index+m)]
        ys = Y[1:(index+m)]
    elseif index + m <= length(X)
        xs = X[(index-m):(index+m)]
        ys = Y[(index-m):(index+m)]
    else
        xs = X[(index-m):end]
        ys = Y[(index-m):end]    
    end
    av_x = sum(xs)/length(xs) 
    av_y = sum(ys)/length(ys)
    
    a = sum((xs .- av_x) .* (ys .- av_y) )/ sum( (xs .- av_x).^2)
    b = av_y - a * av_x
    return a*X[index]+b
end

"""Funkcja, która za pomocą loess znajduje przybliżenie funkcji.
     Uzywamy do odszumiania dzwieków"""
function denoising(Y,m= 10)
    X = 1:length(Y)
    y = [loess(index, X,Y,m) for index in 1:length(X)]
    return y
end

end