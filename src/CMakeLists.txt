ADD_EXECUTABLE( time_fft
                time_fft.cpp
                )
TARGET_LINK_LIBRARIES(time_fft fftw3f minc_io)

ADD_EXECUTABLE( time_correlate
                time_correlate.cpp
                )

TARGET_LINK_LIBRARIES(time_correlate fftw3f minc_io)

INSTALL(TARGETS 
      time_fft
      time_correlate
      RUNTIME DESTINATION bin )

INSTALL(PROGRAMS 
    time_correlate.pl
    temporal_filter.pl
    DESTINATION bin )
