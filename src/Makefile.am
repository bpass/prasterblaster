SUBDIRS = gctp_cpp gtest


noinst_LIBRARIES = librasterblaster.a
librasterblaster_a_SOURCES = projectedraster.cpp \
	         reprojector.cpp \
		 resampler.cpp \
                 rasterchunk.cpp \
                 quadtree.cpp \
		 driver.cpp 
#librasterblaster_a_LDFLAGS = $(all_libraries)  -Lgctp_cpp
librasterblaster_a_LIBADD =  gctp_cpp/libgctp_cpp.a

bin_PROGRAMS = prasterblaster
prasterblaster_CXXFLAGS = -g 
prasterblaster_LDFLAGS = $(all_libraries)  -Lgctp_cpp
prasterblaster_LDADD = librasterblaster.a gctp_cpp/libgctp_cpp.a -lgdal
prasterblaster_SOURCES =  main.cpp 

AM_CPPFLAGS = $(all_includes) -IpRPL -Igctp_cpp

AM_CXXFLAGS = -W -Wall -Wextra -Wcast-align \
		-Wpointer-arith -Wsign-compare -g\
		-Wformat=2 -Wno-format-y2k  -Wmissing-braces -Wparentheses \
		-Wtrigraphs -Wstrict-aliasing=2
