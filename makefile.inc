-include fort.dep

%: %.o
	${F90} ${FFLAGS} ${LDFLAGS} -o ${EXEDIR}/$@ $^

%.o: %.F90
	${F90} ${FFLAGS} -c $<

%.o: %.f90
	${F90} ${FFLAGS} -c $<
