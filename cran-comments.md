# Resubmission

This is a resubmission. 
Thank you for your feedback, here you can find some comments on how I solved the pacakge issues. 

> You do not need to single quote your whole title. Please remove the quotes.

I removed the quotes from the title.

> Please rather use the Authors@R field and declare Maintainer, Authors
and Contributors with their appropriate roles with person() calls.
e.g. something like:
Authors@R: c(person("Alice", "Developer", role = c("aut", "cre","cph"),
                      email = "alice.developer@some.domain.net"),
                     person("Bob", "Dev", role = "aut") )

I changed the DESCRIPTION file to account for authors with this format.

> I think there are too many blank spaces between some words in your
description text. Please check that.

I corrected the description formatting.

> If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
auto-linking.
(If you want to add a title as well please put it in quotes: "Title")

I added the reference on the works from which our package implementation was inspired. 
Unfortunately, there is no reference for the extensions.

> Please add \value to .Rd files regarding exported methods and explain
the functions results in the documentation. Please write about the
structure of the output (class) and also what the output means. (If a
function does not return a value, please document that too, e.g.
\value{No return value, called for side effects} or similar)
Missing Rd-tags:
      clmplus.default.Rd: \value
      clmplus.Rd: \value
      plot.clmplusmodel.Rd: \value
      plot.RtTriangle.Rd: \value
      plotresiduals.clmplusmodel.Rd: \value
      plotresiduals.default.Rd: \value
      plotresiduals.Rd: \value

I added the value tag and output descriptions as raccomended.

# R cmd check 

I obtained 0 ERRORS, 0 WARNINGS and 0 NOTES.