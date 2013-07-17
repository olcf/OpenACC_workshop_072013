default:
	@echo "make clean or make tgz"

clean:
	@cd Practical0; $(RM) myjoblog; cd ..
	@cd Practical1; make clean; cd ..
	@cd Practical2/F; make clean; cd ../..
	@cd Practical2/C; make clean; cd ../..
	@cd Practical3/F; make veryclean; cd ../..
	@cd Practical3/C; make veryclean; cd ../..

tgz: clean
	@cd ..; tar zcf Cray_OpenACC_training_$(shell date +%y%m%d_%H%M%S).tgz --exclude=.git --exclude=TODO --exclude=Makefile_versions --exclude=himeno_?_master.* --exclude=pyver.py --exclude=Advanced --exclude=README_developer.txt --exclude=test*.c Cray_OpenACC_training; cd -

tgz_messmer: clean
	@cd ..; tar zcf Cray_OpenACC_training_$(shell date +%y%m%d_%H%M%S).tgz --exclude=.git --exclude=TODO --exclude=Makefile_versions --exclude=himeno_?_master.* --exclude=pyver.py --exclude=Practical3/C --exclude=Advanced --exclude=README_developer.txt --exclude=test*.c --exclude=Practical1 --exclude=Practical2 --exclude=Practical3 Cray_OpenACC_training; cd -
