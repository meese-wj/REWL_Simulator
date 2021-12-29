add_subdirectory(Square_Periodic_Lattices)

add_library(AddressBooks
            address_books.hpp
            lattice_address_book.hpp)   

set_target_properties(AddressBooks PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(AddressBooks address_books_for_square_periodic_lattices)