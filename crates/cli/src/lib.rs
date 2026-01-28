//! MolMason cli crate

pub fn hello() -> &'static str {
    "MolMason cli"
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_hello() {
        assert!(hello().contains("MolMason"));
    }
}
