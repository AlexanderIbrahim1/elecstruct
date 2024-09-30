# NOTES

When implementing `operator+()` for a type:
  - implement `auto operator+=(const Type& other) -> Type` as a member function
    - it should modify the current instance and return *this (of type `Type&`)
  - implement `auto operator+(Type left, const Type& right) -> Type`
    - it will create a copy of `left`
    - you will modify it with `+=`
