#!/bin/sh

## Note: one can also use ``pdfunite`` instead of ``ghostscript``

#
pdflatex mypart_syllabus_2244_CSCI4650-E01_2024summer
#
gs -dPrinted=false -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=syllabus_2244_CSCI4650-E01_2024summer.pdf -dBATCH mypart_syllabus_2244_CSCI4650-E01_2024summer.pdf 2244_4650_numerical_analysis___course_topics.pdf syllabus_policies_insert.pdf student_services_and_calendar.pdf AcademicCalendar_Summer2024.pdf
#


#
pdflatex mypart_syllabus_2244_MATH4650-E01_2024summer
#
gs -dPrinted=false -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=syllabus_2244_MATH4650-E01_2024summer.pdf -dBATCH mypart_syllabus_2244_MATH4650-E01_2024summer.pdf 2244_4650_numerical_analysis___course_topics.pdf syllabus_policies_insert.pdf student_services_and_calendar.pdf AcademicCalendar_Summer2024.pdf
#

#
pdflatex mypart_syllabus_2244_MATH5660-E01_2024summer
#
gs -dPrinted=false -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=syllabus_2244_MATH5660-E01_2024summer.pdf -dBATCH mypart_syllabus_2244_MATH5660-E01_2024summer.pdf 2244_4650_numerical_analysis___course_topics.pdf syllabus_policies_insert.pdf student_services_and_calendar.pdf AcademicCalendar_Summer2024.pdf
#

return 0

#
# remove the Maymester calendar and 4-week class calendar 
# only needs to be done once
#
gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dFirstPage=1 -dLastPage=2 -sOutputFile=tmp.pdf summer-2022-printable.pdf
mv tmp.pdf summer-2022-printable.pdf

