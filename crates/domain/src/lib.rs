//! MolMason domain crate

pub fn hello() -> &'static str {
    "MolMason domain"
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_hello() {
        assert!(hello().contains("MolMason"));
    }
}
